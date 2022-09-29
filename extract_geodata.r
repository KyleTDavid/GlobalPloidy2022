#script used to extract environmental data for every filtered coordinate, then average them across species 
.libPaths("~/.rlib")

library(dplyr)
library(raster)

#load data
df <- read.csv('Amphibia_coords.txt', sep='\t')

# read rasters
glaciation <- raster("geodata/lgm_global_zero.grd")

temp <- raster("geodata/wc2.1_5m_meantemp.tif")
prec <- raster("geodata/wc2.1_5m_prec.tif")
temp_season <- raster("geodata/wc2.1_5m_seasontemp.tif")
prec_season <- raster("geodata/wc2.1_5m_seasonprec.tif")
alt <- raster("geodata/wc2.1_5m_elev.tif")

delta_temp <- temp - (raster("geodata//cclgmbimeantemp.tif")/10)
delta_prec <- prec - raster("geodata//cclgmbiprec.tif")
delta_temp_season <- temp_season - raster("geodata/cclgmbiseasontemp.tif")
delta_prec_season <- prec_season - raster("geodata/cclgmbiseasonprec.tif")

human_impacts <- raster("geodata/human_impacts/hii_v2geo/dblbnd.adf")
land_use <- raster("geodata/gl-croplands-geotif/cropland.tif") + raster ("geodata/gl-pastures-geotif/pasture.tif")

total_richness <- raster("geodata/species_richness.tif")

points <- SpatialPoints(data.frame(lon=df$longitude, lat=df$latitude))

df$absolute_latitude <- abs(df$latitude)
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

df_mean <- df %>%
  group_by(species) %>%
  summarise_at(vars(3:16), list(mean), na.rm=T)
  
write.table(df_mean, file = "Amphibia_geodata.txt", quote = FALSE, sep = "\t", row.names = FALSE)
