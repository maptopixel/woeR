#' Calc Weights Function
#'
#' This function calculates WoE weights of a raster. Currently doesn't handle missing data.
#' @param evidenceRaster Raster of of evidence
#' @param predictTPts Training points
#' @keywords WoE, Weights of Evidence
#' @export
#' @examples
#' calcWeights()
calcWeights = function (evidenceRaster,predictTPts) {
  r = raster(evidenceRaster)
  #Compute the area for each class
  freqOfRDf = data.frame(freq(r))
  thisNumAreaUnits= freqOfRDf$count

  #get total area
  TotalArea = sum(freqOfRDf$count)

  #Get the values at training pts
  predTpValues = extract(r, predictTPts)

  #summarize the counts using a tabulation
  summaryTableDF = data.frame(table(predTpValues))

  thisNumTP = summaryTableDF$Freq
  classLabels =summaryTableDF$predTpValues

  #total number of TP
  NumTP = sum(summaryTableDF$Freq)

  #classLabels = unique(predTpValues)


  #thisNumTP = c(72,2,2,10,8,6)
  #thisNumAreaUnits = c(5598238, 160302, 126383, 110037, 99895, 79145)


  TotalAreaVector = rep(TotalArea,length(thisNumTP))

  #P(D)
  PD = NumTP/TotalAreaVector

  #P(!D)
  PNotD = (TotalAreaVector - NumTP)/TotalAreaVector

  #Conditional probabilities
  #####################
  #P(B|D)
  PB_D = thisNumTP/NumTP

  #P(B|!D)
  PB_NotD = (thisNumAreaUnits-thisNumTP) / (TotalAreaVector-NumTP)

  #P(!B|D)
  PNotB_D = (NumTP-thisNumTP)/NumTP

  #P(!B|!D)
  PNotB_NotD = (TotalAreaVector-thisNumAreaUnits-NumTP-thisNumTP)/(TotalAreaVector-NumTP)

  #Calculate the weights
  #####################
  #Wplus = LN(P(B|D) / P(B|!D))
  Wplus = log(PB_D/PB_NotD)

  #Wminus = LN(P(!B|D) / P(!B|!D))
  Wminus = log(PNotB_D/PNotB_NotD)


  #Determine new class values
  newClassValues = Wplus #This needs some checking based on a cumululative asceding / descending param

  revisedClassesDf = data.frame(cbind(classLabels,newClassValues))
  names(revisedClassesDf) = c("id","v")

  #subsitute in the weight values into the class values
  weightedClassRaster <- subs(r, revisedClassesDf)

  return(weightedClassRaster)
}
