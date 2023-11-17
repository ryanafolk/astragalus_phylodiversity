# astragalus_phylodiversity
## Base-level directory
`astragalus.R` is the main analysis script. The base-level directory also provides two map summaries at differing spatial grain: `Astragalus_biodiverse_update_15km.pdf` shows 15-km square grid cells, and `Astragalus_biodiverse_update_50km.pdf` shows the same at 50-km square grid cell size, with the latter used for downstream analysis.

## Astragalus_50km_tifsToShare and astragalus_mikebelitz_15km_outputs
These two directories give a series of diversity statics globally at the two spatial grains (respectively in each folder) in GeoTIFF format. For instance, `astragalus_PHYLO_RPD2.tif` represents RPD.

## Astragalus_50km_csvsToShare and Astragalus_CSVs_ToShare_15km
These two directories represent the same data as above, but in tabular CSV format: x, y, value.

## environmental_data
This directory gives environmental data extracted for grid cells, using the 50 km grid cells: x, y, value.

## R_figs
This directory gives a series of plots of various environmental parameters against the diversity metrics. File names indicate the comparison with choices of plotting following the narative of the main text. For instance, `rpd_latitude_NORTHAMERICA.pdf` shows RPD vs. latitude, showing only grid cells for North America.

## new_figs
This directory gives prepared map figures.

