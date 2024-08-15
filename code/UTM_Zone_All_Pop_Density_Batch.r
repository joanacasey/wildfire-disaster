# Description -------------------------------------------------------------

# Author: Alex Northrop, Benjamin Steiger
# Author: Joan Casey, Lauren Wilner, and Brian High 
# Date: 10/31/2022
# Last Update: 07/15/2024
# Goal: Population Density Overlap Batch

# Load packages -----------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(raster, tmap, SpatialKDE, dplyr, here, sf, readr, beepr, 
               parallel, stringr, purrr, furrr)


# Setup the parallel processing plan
#plan(multisession(workers = 64))

# Load functions ----------------------------------------------------------
pop_den_work <- function(.data, pop_raster, crs) {
  fire_work <- .data %>% 
    st_transform(crs = proj4string(pop_raster))
  
  disaster_id_obj <- fire_work$disaster_id
  
  # Crop and mask for wildfires 
  
  r <- crop(pop_raster, extent(fire_work))
  r <- mask(r, fire_work)
  
  ## Convert Raster to Points 
  
  utm_10_points = rasterToPoints(r, fun=function(x){x>=1}, spatial = TRUE)
  
  utm_10_sf <- utm_10_points %>% 
    st_as_sf(coords = c("x", "y"), dim = "XY") 
  
  # transform to UTM zone NAD83 meters unit
  
  utm_10_sf$geometry <- st_transform(utm_10_sf$geometry, crs = crs)
  fire_work$shape <- st_transform(fire_work$shape, crs = crs)
  
  cell_size <- 300
  band_width <- 907 # Radius in meters of a mile^2 
  
  test <- as.data.frame(utm_10_sf)
  
  if(nrow(test) > 0) {
    
    #Check if this keeps just points that joint later
    #save disaster id
    fire_name_utm_10 <- utm_10_sf %>% 
      st_join(fire_work)%>%
      dplyr::select(disaster_id)

    #create grid
    grid_utm_10 <- fire_name_utm_10 %>% 
      create_grid_rectangular(cell_size = cell_size, side_offset = band_width) 
    
    #keep just one row disaster id
    fire_name_utm_10 <- fire_name_utm_10[1,1] 
    
    density <- utm_10_sf %>% 
      kde(
        band_width = band_width, 
        kernel = "quartic", 
        grid = grid_utm_10, 
        weights = utm_10_sf[[1]]
    )
    
    density_criteria <- as.data.frame(density) %>% 
      filter(kde_value == max(kde_value)) %>% 
      dplyr::distinct(kde_value) %>%
      dplyr::mutate(density_criteria_met = ifelse(kde_value >= 250, TRUE, FALSE))
    
    density_criteria$disaster_id <- fire_name_utm_10
    density_criteria
    
  }
  # else{
  #   
  #   density_criteria <- utm_10_sf %>% mutate(kde_value=0, density_criteria_met = FALSE, disaster_id = disaster_id_obj)  
  #   
  # }
  
  #return(density_criteria)
  
  
}

pop_den_check <- function(df.in, pop_raster, crs = crs_espg) {
    future_map_dfr(.x = 1:nrow(df.in), ~ pop_den_work(.data = df.in[.x, ], 
                                                      pop_raster = pop_raster, crs = crs)) 

    #print(i)
}

# Load data ---------------------------------------------------------------

# binded all disaster fires, CONUS, 2000-2019

load(here(
  "data",
  "raw", 
  "all_disasters_select_vars.rdata"
))


# clarify these datasets
conus <- all_disaster_perimeters_buffers_conus_dist_select_vars
hawaii <- all_disaster_perimeters_buffers_hawaii_dist_select_vars
alaska <- all_disaster_perimeters_buffers_alaska_dist_select_vars

# fix the disaster id variable
conus <- conus %>% mutate(disaster_id = sub(";.*", "", disaster_nested_id))
hawaii <- hawaii %>% mutate(disaster_id = sub(";.*", "", disaster_nested_id)) #in espg=2784
alaska <- alaska %>% mutate(disaster_id = sub(";.*", "", disaster_nested_id)) #in espg=3338

# acerage
conus <- conus %>% mutate(st_area_acre_final = area_sq_m*.000247105)
hawaii <- hawaii %>% mutate(st_area_acre_final = area_sq_m*.000247105)
alaska <- alaska %>% mutate(st_area_acre_final = area_sq_m*.000247105)

# crs UTM mapping by place
utm_crs <- read.csv("data/utm_popden.csv")

# Turn off spherical geometry ---------------------------------------------

sf_use_s2(FALSE)

# Subset fires to those without empty geometries -------------------------
# empty geometries won't work

# rename df
#fires <- binded_conus_disaster_fires_2000_2019
lst <- lapply(list(conus, hawaii, alaska), function(fires) {
  #fires <- alaska[1:2,]
  
  # subset to non-empty geometries
  fires <- fires %>%
    filter(!is.na(perimeter_source))
  
  # Subset to appropriate UTM  --------------------------------------------
  
  # UTM zone
  # State list
  lst.utm <- lapply(1:nrow(utm_crs), function(i) {
    utm_crs_i <- utm_crs[i,]
    
    u <- utm_crs_i$utm

    utm_crs_10_temp <- utm_crs_i %>%
      mutate(states_list = str_split(states, ", "))
    
    states_vector <- unlist(utm_crs_10_temp$states_list[1])
    
    fires_utm_10 <- fires #%>%
    #  filter(disaster_list_states %in% states_vector)
    # debug 
    
    if(nrow(fires_utm_10)>0){
    
      crs_espg <- utm_crs_10_temp$crs
      # print(unique(crs_espg))
      # print(unique(fires_utm_10$disaster_list_states))
      
      # rasters
      us_raster_00 <- raster(file.path("data/casey-cohort", utm_crs_i$pop_raster_00), crs = st_crs(4326)$wkt)
      us_raster_10 <- raster(file.path("data/casey-cohort", utm_crs_i$pop_raster_10), crs = st_crs(4326)$wkt)
      
      
      # Create the buffer for wildfire boundaries -------------------------------
      
      # # change crs to us_raster_00, WGS84 for HI and AK
      # 
      # fires_utm_10 <- fires_utm_10 %>%
      #   st_transform(crs = proj4string(us_raster_00))
      
      # add sizing categories - large = > 1000 acres, small = ≤ 1000
      # buffer on new projected shape
      
      fires_utm_10 <- fires_utm_10 %>%
        mutate(
          size = ifelse(st_area_acre_final > 1000, "large_fire", "small_fire"),
          shape = st_buffer(shape, ifelse(size == "large_fire", 20000, 10000))
        )
      
      # change shape to character, select shape and disaster_id
      
      fire_geom_id <- as.data.frame(fires_utm_10) %>%
        mutate(shape = as.character(shape)) %>%
        dplyr::select(disaster_id, shape)
      
      # subset to 2000-2009 fires
      
      fires_utm_10_00 <- fires_utm_10 %>%
        filter(as.numeric(aggregate_year) <= "2009")
      
      # subset to 2010-2019 fires_utm_10
      
      fires_utm_10_10 <- fires_utm_10 %>%
        filter(as.numeric(aggregate_year) >= "2010")
      
      # create 2000 and 2010 empty output dataframes
      
      #output_00 <- data.frame()
      #output_10 <- data.frame()
      
      ### 2000-2009 Fires
      
      
      # 2000-2009 fires ---------------------------------------------------------
      # df.in <- fires_utm_10_00[1:2,]
      # pop_raster <- us_raster_00
      # output_00 <- pop_den_check(test_df, output_00, us_raster_00)
      # .data = fires_utm_10_00
      # pop_raster = us_raster_00
      # crs = 26912
      
      output_00 <- pop_den_check(fires_utm_10_00, us_raster_00, crs=crs_espg)
      output_10 <- pop_den_check(fires_utm_10_10, us_raster_10, crs=crs_espg)
      
      # Output 2000-2009 --------------------------------------------------------
      
      # create 2000-09 output dataframe, filter out duplicate fires by disaster_id, drop geometry
      output_00 <- as.data.frame(output_00) %>%
        filter(!duplicated(disaster_id)) %>%
        dplyr::select(-geometry)
      
      # Output 2010-19 fires ----------------------------------------------------
      
      # create 2010-19 output dataframe, filter out duplicate fires by disaster_id, drop geometry
      output_10 <- as.data.frame(output_10) %>%
        filter(!duplicated(disaster_id)) %>%
        dplyr::select(-geometry)
      
      # Total output ------------------------------------------------------------
      
      # bind 2000-09 and 2010-19 outputs
      rbind(output_00, output_10)
      
      
        
    }
  })
  rbind(lst.utm)
})

tot_output <- rbind(lst)

# write total output density criteria .csv
tot_output %>%
  write_csv(
    here(
      "data",
      "processed",
      "combined",
      "pop_density_test",
      "density_criteria_all.csv"
    )
  )
