#! /usr/bin/env bash

mkdir -p $(dirname "$0")/../data/raw/pop_data
cd $(dirname "$0")/../data/raw/pop_data 

# Download raw zips
wget -v https://socscape.edu.pl/socscape_data/us_grids/us_pop/us_pop2000myc.zip
wget -v https://socscape.edu.pl/socscape_data/us_grids/us_pop/us_pop2010myc.zip
wget -v https://socscape.edu.pl/socscape_data/us_grids/us_pop/us_pop2020myc.zip
wget -v https://data.worldpop.org/GIS/Population/Global_2000_2020_1km/2000/USA/50_US_states_1km_2000.zip
wget -v https://data.worldpop.org/GIS/Population/Global_2000_2020_1km/2010/USA/50_US_states_1km_2010.zip
wget -v https://data.worldpop.org/GIS/Population/Global_2000_2020_1km/2020/USA/50_US_states_1km_2020.zip

# Extract full 30m zips
unzip us_pop2000myc.zip
unzip us_pop2010myc.zip
unzip us_pop2020myc.zip

# Extract needed files from 1km zips
unzip -j 50_US_states_1km_2000.zip US-HI_ppp_2000_1km.tif  
unzip -j 50_US_states_1km_2000.zip US-AK_ppp_2000_1km.tif 
unzip -j 50_US_states_1km_2010.zip US-HI_ppp_2010_1km.tif  
unzip -j 50_US_states_1km_2010.zip US-AK_ppp_2010_1km.tif 
unzip -j 50_US_states_1km_2020.zip US-HI_ppp_2020_1km.tif  
unzip -j 50_US_states_1km_2020.zip US-AK_ppp_2020_1km.tif 

# Remove unneeded zips
rm -f us_pop20*myc.zip
rm -f 50_US_states_1km_20*.zip 
