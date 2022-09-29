#script used to further filter and format GBIF coordinates for downstream analyses 
library(rotl)
library(dplyr)
library(tidyr)
library(ape)
library(countrycode)
library(CoordinateCleaner)

args = commandArgs(trailingOnly=TRUE)

#input is GBIF coordinates file, for example https://doi.org/10.15468/dl.2jxpma
df <- read.csv(args[1], sep='\t', header=FALSE, comment.char='', quote="")
colnames(df) <- c("gbifID", "datasetKey", "occurrenceID", "kingdom", "phylum", "class", "order", "family", "genus", "species", "infraspecificEpithet", "taxonRank", "scientificName", "verbatimScientificName", "verbatimScientificNameAuthorship", "countryCode", "locality", "stateProvince", "occurrenceStatus", "individualCount", "publishingOrgKey", "decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters", "coordinatePrecision", "elevation", "elevationAccuracy", "depth", "depthAccuracy", "eventDate", "day", "month", "year", "taxonKey", "speciesKey", "basisOfRecord", "institutionCode", "collectionCode", "catalogNumber", "recordNumber", "identifiedBy", "dateIdentified", "license", "rightsHolder", "recordedBy", "typeStatus", "establishmentMeans", "lastInterpreted", "mediaType", "issue")

# if species name is missing - take it from "scientific name"
ind = which(is.na(df$species) & !is.na(df$verbatimScientificName))
df$species[ind] = df$verbatimScientificName[ind]

# remove rows without species name
df = df[which(!is.na(df$species)),] 

# filter out unnecessary columns 
df <- df %>% select(species, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, countryCode)

# resolve names with OTT
resolved_names <- tnrs_match_names(unique(df$species), "Amphibians")
df$species <- resolved_names[match(tolower(df$species), resolved_names$search_string) , 'unique_name']
df <- df %>% drop_na(species, decimalLatitude, decimalLongitude)
df$species <- gsub(" \\(species in.*", "", df$species)
#trim subspecies
df$species <- sub("^([^ ]* [^ ]*).*", "\\1", df$species)

# filter entries not in tree
tree <- read.tree('~/Desktop/global_poly/Amphibia_tree.nh')
df <- df[gsub(" ", "_", df$species) %in% tree$tip.label,]

# convert country code from ISO2c to ISO3c
df$countryCode <- countrycode(df$countryCode, origin = 'iso2c', destination = 'iso3c')
# convert latitude, longitude and coordinateUncertainty to numeric values
df[,c(2:4)] = apply(df[,c(2:4)], 2, function(x) as.numeric(as.character(x)))

# filter entries if reported coordinate uncertainty is greater than 100km
df <- df[is.na(df$coordinateUncertaintyInMeters) | df$coordinateUncertaintyInMeters < 100000,]

# resolution of coordinates must be at least 0.1
length_decimal_lat = nchar(sub("^[^.]*", "", df$decimalLatitude))-1
length_decimal_lon = nchar(sub("^[^.]*", "", df$decimalLongitude))-1
df = df[which(length_decimal_lat>=1 & length_decimal_lon>=1),] 

# clean coordinates
flags <- clean_coordinates(x = df,
                    lon = "decimalLongitude",
                    lat = "decimalLatitude",
                    countries = "countryCode",
                    species = "species",
                    tests = c("capitals", "centroids", "equal", "gbif", "institutions", "zeros")
)
df <- df[flags$.summary,] 
df <- df %>% select(species, decimalLatitude, decimalLongitude)

#round to nearest 0.1 decimal and take unique value
df$decimalLatitude <- round(df$decimalLatitude ,digit=1)
df$decimalLongitude <- round(df$decimalLongitude ,digit=1)

df <- unique(df)

df <- rename(df, latitude=decimalLatitude, longitude=decimalLongitude)
df <- na.omit(df)
write.table(df, file = "Amphibia_coords.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

