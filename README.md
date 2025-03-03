# herbarium_project

Collection of scripts and notes for the herbarium projects.


# Data collection
For every occurrence of Amaranthus tuberculatus with coordinates in the GBIF data, we reprojected the coordinates to the coordinate reference system (CRS) of the [Annual National Land Cover Database (NLCD)](https://www.usgs.gov/centers/eros/science/annual-national-land-cover-database). Then, based on the year the specimen was collected, we referenced that year's NLCD data to extract the land class value of the exact pixel the coordinates corresponded to.
For occurrences in the ambiguous "Barren Land", we cross-referenced them with Google Maps Street View, aerial imagery, and 3D imagery to manually classify whether they were in human-mediated or natural areas.
