# Description -------------------------------------------------------------

# Author: Alex Northrop, Benjamin Steiger
# Date: 10/31/2022
# Last Update: 08/25/2023
# Goal: Population Density Overlap Batch - UTM Zone 10

# Load packages -----------------------------------------------------------
library(pacman)
p_load(raster, tmap, SpatialKDE, dplyr, here, sf, readr, beepr, parallel, stringr)
# test

# Load functions ----------------------------------------------------------
pop_den_check <- function(df.in, df.out, pop_raster, crs = crs_espg) {

    lst <- mclapply(1:nrow(df.in), function(i) {
    fire_work <- df.in[i,] %>% 
      st_transform(crs = proj4string(pop_raster))
    disaster_id_obj <- fire_work$disaster_id
    
    # Crop and mask for UTM Zone 10 wildfires 
    
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
    
  if(sum(test$us_pop2000myc.population_data) > 0) {
    
      grid_utm_10 <- utm_10_sf %>% 
        st_join(fire_work) %>% 
        create_grid_rectangular(cell_size = cell_size, side_offset = band_width) 
      
      density <- utm_10_sf %>% 
        kde(band_width = band_width, kernel = "quartic", grid = grid_utm_10, weights = utm_10_sf$us_pop2000myc.population_data)
      
      density_criteria <- as.data.frame(density) %>% 
        filter(kde_value == max(kde_value)) %>% 
        dplyr::mutate(density_criteria_met = ifelse(kde_value >= 250, "Yes","No"))
      
      density_criteria$disaster_id <- fire_work$disaster_id
    
      
  } else{
      density_criteria <- utm_10_sf %>% mutate(kde_value=0, density_criteria_met = "No", disaster_id = disaster_id_obj)  
      
      }
    }, mc.cores = 4 )
    
    df.out <- bind_rows(df.out, lst)
    
    #print(i)
}

# Load data ---------------------------------------------------------------

# binded all disaster fires, CONUS, 2000-2019

load(here(
  "data",
  "raw",
  "all_disasters_select_vars.rdata"
))

load(here(
  "data",
  "processed",
  "combined",
  "binded_conus_disaster_fires_2000_2019.rdata"
))

# clarify these datasets
conus <- all_disaster_perimeters_buffers_conus_dist_select_vars
hawaii <- all_disaster_perimeters_buffers_hawaii_select_vars
alaska <- all_disaster_perimeters_buffers_alaska_dist_select_vars

# fix the disaster id variable
conus <- conus %>% mutate(disaster_id = sub(";.*", "", disaster_nested_id))
hawaii <- hawaii %>% mutate(disaster_id = sub(";.*", "", disaster_nested_id))
alaska <- alaska %>% mutate(disaster_id = sub(";.*", "", disaster_nested_id))

# rasters
us_raster_00 <- raster("~/casey-cohort/us_pop2000myc.tif")
us_raster_10 <- raster("~/casey-cohort/us_pop2010myc.tif")
# us_raster_00 <- raster("/Volumes/casey-cohort/us_pop2000myc.tif")
# us_raster_10 <- raster("/Volumes/casey-cohort/us_pop2010myc.tif")

# crs UTM mapping by place
utm_crs <- read.csv(here("data", "raw", "utm_popden.csv"))

# Turn off spherical geometry ---------------------------------------------

sf_use_s2(FALSE)

# Subset fires to those without empty geometries -------------------------

# empty geometries won't work

# rename df
fires <- binded_conus_disaster_fires_2000_2019

# subset to non-empty geometries
fires <- fires %>%
  filter(!is.na(perimeter_source))

# Subset to appropriate UTM  --------------------------------------------

# UTM zone
# State list
for(u in utm_crs$utm){
  
  print(u)

  utm_crs_10_temp <- utm_crs %>%
    filter(utm==u) %>%
    mutate(states_list = str_split(states, ", "))
  states_vector <- unlist(utm_crs_10_temp$states_list[1])

  fires_utm_10 <- fires %>%
    filter(state %in% states_vector)
  
  crs_espg <- utm_crs_10_temp$crs
  print(unique(crs_espg))
  print(unique(fires_utm_10$state))

# Create the buffer for wildfire boundaries -------------------------------

# change crs to us_raster_00

fires_utm_10 <- fires_utm_10 %>%
  st_transform(crs = proj4string(us_raster_00)) 

# add sizing categories - large = > 1000 acres, small = ≤ 1000
# buffer on new projected shape

fires_utm_10 <- fires_utm_10 %>%
  mutate(size = ifelse(st_area_acre_final > 1000, "large_fire", "small_fire"),
         shape = st_buffer(shape, ifelse(size == "large_fire", 20000, 10000)))

# change shape to character, select shape and disaster_id

fire_geom_id <- as.data.frame(fires_utm_10) %>%
  mutate(shape = as.character(shape)) %>%
  dplyr::select(disaster_id, shape)

# subset to 2000-2009 fires

fires_utm_10_00 <- fires_utm_10 %>%
  filter(as.numeric(year) <= "2009")

# subset to 2010-2019 fires_utm_10

fires_utm_10_10 <- fires_utm_10 %>%
  filter(as.numeric(year) >= "2010")

# create 2000 and 2010 empty output dataframes

output_00 <- data.frame()
output_10 <- data.frame()

### 2000-2009 Fires


# 2000-2009 fires ---------------------------------------------------------
# df.in <- fires_utm_10_00[1:2,]
# pop_raster <- us_raster_00
# output_00 <- pop_den_check(test_df, output_00, us_raster_00)

output_00 <- pop_den_check(fires_utm_10_00, output_00, us_raster_00)
output_10 <- pop_den_check(fires_utm_10_10, output_10, us_raster_10)

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
tot_output <- rbind(output_00, output_10)

# write total output density criteria .csv
tot_output %>%
  write_csv(here("data", "processed", "combined", "pop_density_test", "utm_zone_10", "density_criteria.csv"))

#join fires to total output by disaster_id
fire_density_criteria_full <- fires_utm_10_00 %>%
  left_join(tot_output, by = "disaster_id")

# write full fire density criteria .csv file
fire_density_criteria_full %>%
  write_csv(here("data", "processed", "combined", "pop_density_test", paste0("utm_zone_", utm), "fire_density_criteria_full.csv"))
}
