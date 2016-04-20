# woer R package #
### Some R functions for computing weights of evidence##

Weights of evidence caluclation for spatial data following Bonham-Carter, 1994.


* Determine the weights of evidence for a classified raster using a file (e.g. shp) of training points:

`predictTPts<- readOGR(dsn="your_dir", layer="training_points_sample")`

`evidenceRasterName= paste0("your_dir","/", "classified_raster.tif"); `

`isCumulative = TRUE`

`evidenceRaster = readGDAL(evidenceRasterName)`

`weightedEvidenceRaster= calcWeights(evidenceRaster,predictTPts,isCumulative)`



