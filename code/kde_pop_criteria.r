if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(raster, tmap, SpatialKDE, dplyr, here, sf, readr, parallel, stringr, purrr, furrr, terra, tidyr, logger, glue, fs)

## Package Options
options("readr.show_col_types" = FALSE)
rasterOptions(
    maxmemory = 5e09, 
    chunksize = 1e09, 
    todisk = FALSE,
    tmpdir = dir_create('~/rtmp', recurse = TRUE) # Use a non-default spot on brain for temp file storage -- the default has less space than home directory
)
log_appender(appender_tee('fire.log')) # log to file

## Future Plan 

#plan( # Switch on to use batchtools for Sun Grid Engine (will have to tweak for Hyak)
#    tweak(batchtools_sge, template = 'sge.tmpl', resources = list(ncpus = 8))
#)

plan(sequential) # Switch on to run sequentially for testing

# plan(multisession, workers = 48) # switch on to run on single node

## Logging and Saving Partial Results 

# Partial results saved under PID of main process
runid <- Sys.getpid()
dir.create(glue('data/partial/{runid}/'), recursive = TRUE)
log_info('Begin Run - PID: {runid}')

## Process Data

load('data/raw/all_disasters_select_vars.rdata')

# Read in fires and adjust some columns 
log_info('Reading in fire perimeters.')
fires <- bind_rows(
    hawaii = st_transform(all_disaster_perimeters_buffers_hawaii_dist_select_vars, 4326)#,
    #alaska = st_transform(all_disaster_perimeters_buffers_alaska_dist_select_vars, 4326),
    #conus = st_transform(all_disaster_perimeters_buffers_conus_dist_select_vars, 4326)
) %>% 
    mutate(disaster_id = sub(";.*", "", disaster_nested_id)) %>% 
    mutate(st_area_acre_final = area_sq_m*.000247105) %>% 
    mutate(size = ifelse(st_area_acre_final > 1000, "large_fire", "small_fire")) 
    
# Assign working CRS and decade, split by that plus the buffer that will be used
utm_crs <- read_csv("data/utm_popden.csv") %>%
    pivot_longer(matches('pop_raster'), names_to = 'pop_raster_year', values_to = 'pop_raster_file') %>%
    mutate(pop_raster_file = file.path('data/raw/pop_data', pop_raster_file)) %>%
    mutate(pop_raster_year = str_extract(pop_raster_year, '[0-9]{2}'))
fires <- fires %>% 
    cross_join(utm_crs) %>% 
    filter(str_detect(states, str_extract(aggregate_states, '^[A-Z]{2}'))) %>% # Match state to work CRS. If fire crosses states, take the first. 
    filter((aggregate_year < 2010 & pop_raster_year == '00') | (aggregate_year >= 2010 & pop_raster_year == '10')) %>% # match year to pop raster
    filter(!is.na(perimeter_source)) %>% # eliminate empty geoms
    mutate() %>% 
    group_by(
        crs, 
        pop_raster_file, 
        buffer_dist = ifelse(size == "large_fire", 20000, 10000)
    ) %>% 
    nest()

# Project each group into its working CRS, buffer the fire geom
log_info('Projecting fire perimeters and calculating buffers.')
fires <- fires %>% 
    mutate(data = map(data, st_transform, crs)) %>% 
    mutate(data = map(data, st_buffer, buffer_dist)) %>% # apply buffer to fire perimeters
    ungroup()

# Do some preprocessing of rasters so that they aren't full size for every single fire
# !!!!!! This takes way too long outside HI and AK !!!!!!!!!!!!!!!!!
log_info('Preprocessing population rasters.')
pop_rasters <- fires %>% 
    select(pop_raster_file, data) %>%
    future_pmap(
        function(pop_raster_file, data){
            log_appender(appender_file('fire.log'))
            log_info('Preprocessing population rasters for projection: {st_crs(data)$epsg}')    
            raster(pop_raster_file) %>%
                projectRaster(crs = st_crs(data)$wkt) %>%
                crop(extent(st_bbox(data))) %>% 
                mask(data)
        },
        .options = furrr_options(seed = nchar('dummy'))
    )
names(pop_rasters) <- fires$pop_raster_file
write_rds(pop_rasters, glue('data/partial/{runid}/__preprocessed_pop_rasters.rds'))

# Break out all individual fires into a list of dfs (a bit messy since
#  no different CRSes can live in the same sf data frame)
fires <- fires %>% 
    group_split(row_number(), .keep = FALSE) %>% 
    map(unnest, cols = c(data)) %>%
    map(function(crs) group_split(crs, row_number(), .keep = FALSE)) %>%
    unlist(recursive = FALSE) %>%
    map(st_as_sf)

# Add preprocessed fires
fires <- map(
    1:length(fires), 
    function(i){
        list(fire = fires[[i]], pop_raster = pop_rasters[[fires[[i]]$pop_raster_file]])
    }
)

# Loop through individual fires determining criteria for each
log_info('Calculating Population KDEs.')
density_criteria <- future_map(
    fires,
    function(fire){
        tryCatch(
            {   
                fire_perim <- fire[['fire']]
                pop_raster <- fire[['pop_raster']]
                log_appender(appender_file('fire.log'))
                pop_raster <- crop(pop_raster, extent(fire_perim))
                pop_raster <- mask(pop_raster, fire_perim)
                pop_sf <- rasterToPoints(pop_raster, fun=function(x){x>=1}, spatial = TRUE) %>% 
                    st_as_sf(coords = c("x", "y"), dim = "XY")

                band_width <- 907 # Radius in meters of a mile^2 
                cell_size <- 300

                fire_grid <- pop_sf %>%
                    st_join(fire_perim) %>%
                    dplyr::select(disaster_id)

                kde_grid <- create_grid_rectangular(
                    fire_grid,
                    cell_size = cell_size, 
                    side_offset = band_width
                ) 

                density <- kde(
                    pop_sf,
                    band_width = band_width, 
                    kernel = "quartic", 
                    grid = kde_grid,
                    weights = pop_sf[[1]]
                )

                density_criteria <- as.data.frame(density) %>% 
                    filter(kde_value == max(kde_value)) %>% 
                    dplyr::distinct(kde_value) %>%
                    dplyr::mutate(density_criteria_met = ifelse(kde_value >= 250, TRUE, FALSE))
                    
                density_criteria$disaster_id <- fire_perim$disaster_id
                log_info('Fire {fire$disaster_id} succeeded.')
                write_rds(density_criteria, glue('data/partial/{runid}/{fire_perim$disaster_id}.rds'))
                density_criteria
            },
            error = function(e){
                log_appender(appender_file('fire.log'))
                log_error('Fire {fire_perim$disaster_id} failed with error code: {e}')
                return(NULL)
            }
        )
    },
    .options = furrr_options(globals = "runid", packages = c("logger", "raster", "dplyr", "SpatialKDE", "sf", "glue", "terra", "readr"))
)

density_criteria <- bind_rows(density_criteria)

write_csv(density_criteria, 'data/processed/density_criteria.csv')

log_info('End Run')
