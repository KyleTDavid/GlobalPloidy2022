# GlobalPloidy2022

#### Scripts to replicate "Global Gradients in the Distribution of Animal Polyploids"

- `filter_coordinates.r` Perform additional filtering steps on GBIF download file [example](https://doi.org/10.15468/dl.2jxpma)

- `extract_geodata.r` Create [environmental data tables](https://figshare.com/articles/dataset/Environmental_Data_Tables/20457030) from filtered coordinates and geographic raster files

- `primary_analyses.r` Run models and generate summary statistics for all primary analysis, using environmental data tables and [phylogenies](http://doi.org/10.6084/m9.figshare.20457045) 

Note that scripts refer to "Amphibia" but can also be used for the other clades by simply changing file names, identical procedures were used for all three groups
