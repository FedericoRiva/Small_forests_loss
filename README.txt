Riva, Martin, Millard and Fahrig (2022)
Loss of the world's smallest forests
Global Change Biology

Data for "Loss of the world's smallest forests" (Riva et al. 2022, Global Change Biology ) is available in Dryad with DOI https://doi.org/10.5061/dryad.69p8cz956

The repository contains the R script and data necessary to replicate the findings presented in Riva et al. (2022) in Global Change Biology.

Running the R script "Riva_et_al_GCB.R" requires opening the following 11 files included in the repository:

Six raster files, all projected in Equal Earth coordinate system
lc1992_waterna.tif - land cover from the year 1992
lc2020_waterna.tif - land cover from the year 2020
biome_rast30.tif - biomes of the Earth
agg1992_30.tif - forest aggregated at a 30 km resolution for visualization, from the year 1992
agg2020_30.tif - forest aggregated at a 30 km resolution for visualization, from the year 2020
lc2020_humandominated.tif - land cover reclassified into anthropogenic vs natural, from the year 2020

Five text files, including the results of randomized point creation across the Earth and data preparation for analysis
data_forest1992.csv - random points created in forests in 1992
data_forest2020.csv - random points created in forests in 2020 (not used for modeling exercise, just statistics comparing forests in 1992 and 2020)
results1992.csv - extraction of forest patch area for points created on 1992 forest raster
results2020.csv - extraction of forest patch area for points created on 2020 forest raster
results_final.csv - table used for models
