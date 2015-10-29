#' Calc Posterior Map Function
#'
#' This function calculates the posterior map based on weights of evidence rasters. The function returns a raster object.
#' @param evidenceRasters (RasterBrick) are multiple classifier rasters according to weighted evidence
#' @keywords WoE, Weights of Evidence
#' @export 
#' @examples
#' calcWeights()
calcPosteriorMap = function (evidenceRasters) {

	#compute prior logit value
	
	#create raster of prilogit value, priLogitRast equal to evidenceRasters dimension
	
	#posteriorLogitRaster = map calc sum of priLogitRast and evidenceRasters

	#posterOdds = exp(posteriorLogitRaster)
	
	#posteriorProb = posterOdds/(1+posterOdds)
}
