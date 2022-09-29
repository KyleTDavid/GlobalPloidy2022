#script used to extract environmental data for every filtered coordinate, then average them across species 
.libPaths("~/.rlib")

#dependencies
library(dplyr)
library(raster)

#load data
df <- read.csv('Amphibia_coords.txt', sep='\t')

#read WorldClim environmental rasters
glaciation <- raster("geodata/lgm_global_zero.grd")
temp <- raster("geodata/wc2.1_5m_meantemp.tif")
prec <- raster("geodata/wc2.1_5m_prec.tif")
temp_season <- raster("geodata/wc2.1_5m_seasontemp.tif")
prec_season <- raster("geodata/wc2.1_5m_seasonprec.tif")
alt <- raster("geodata/wc2.1_5m_elev.tif")

#subtract environmental variables from Community Climate System Model 4 from contemporary estimates to calculate change from LGM
delta_temp <- temp - (raster("geodata//cclgmbimeantemp.tif")/10)
delta_prec <- prec - raster("geodata//cclgmbiprec.tif")
delta_temp_season <- temp_season - raster("geodata/cclgmbiseasontemp.tif")
delta_prec_season <- prec_season - raster("geodata/cclgmbiseasonprec.tif")

#read SEDAC anthropogenic impact rasters
human_impacts <- raster("geodata/human_impacts/hii_v2geo/dblbnd.adf")
land_use <- raster("geodata/gl-croplands-geotif/cropland.tif") + raster ("geodata/gl-pastures-geotif/pasture.tif")

#read custom raster of species richness (see manuscript for details)
total_richness <- raster("geodata/species_richness.tif")

#convert coordinates to SpatialPoints object 
points <- SpatialPoints(data.frame(lon=df$longitude, lat=df$latitude))

#convert latitude to absolute latitude
df$absolute_latitude <- abs(df$latitude)

#for every coordinate extract data point from each raster
df$glaciation <- extract(glaciation, points)
df$temperature <- extract(temp, points)
df$precipitation <- extract(prec, points)
df$temp_season <- extract(temp_season, points)
df$prec_season <- extract(prec_season, points)
df$altitude <- extract(alt, points)
df$delta_temp <- extract(delta_temp, points)
df$delta_prec <- extract(delta_prec, points)
df$delta_ts <- extract(delta_temp_season, points)
df$delta_ps <- extract(delta_prec_season, points)
df$human_impacts <- extract(human_impacts, points)
df$land_use <- extract(land_use, points)
df$total_richness <- extract(total_richness, points)

#get mean value of each variable for every species averaged across all occurrences
df_mean <- df %>%
  group_by(species) %>%
  summarise_at(vars(3:16), list(mean), na.rm=T)

#write output
write.table(df_mean, file = "Amphibia_geodata.txt", quote = FALSE, sep = "\t", row.names = FALSE)
