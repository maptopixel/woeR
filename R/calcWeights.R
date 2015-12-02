#' Calc Weights Function
#'
#' This function calculates WoE weights of a raster. Currently doesn't handle missing data.
#' @param evidenceRaster Raster of of evidence
#' @param predictTPts Training points
#' @keywords WoE, Weights of Evidence
#' @export
#' @examples
#' calcWeights()
calcWeights = function (evidenceRaster,predictTPts,cumulative) {
  r = raster(evidenceRaster)
  #total number of TP
  NumOfTP = length(predictTPts)

  #Compute the area for each class
  freqOfRDf = data.frame(freq(r))
  classIdsList = freqOfRDf$value
  thisNumAreaUnits = freqOfRDf$count

  #Determine if doing cumulative ascending or categorical
  if (cumulative==TRUE) {
    thisNumAreaUnits= cumsum( thisNumAreaUnits)
  }
  thisNumAreaUnits

  #get total area
  TotalArea = sum(freqOfRDf$count)

  #Get the values at training pts
  predTpValues = extract(r, predictTPts)

  #summarize the counts of tps using a tabulation
  summaryTableDF = data.frame(table(predTpValues))
  summaryTableDF

  #Determine if doing cumulative ascending or categorical
  if (cumulative==TRUE) {
    summaryTableDF$Freq = cumsum(summaryTableDF$Freq)
  }


  df = data.frame(classIdsList, rep(NumOfTP, length(classIdsList)))
  colnames(df) = c("ID","Freq")
  colnames(summaryTableDF) = c("ID","Freq")
  summaryTableDFFinal = merge(df,summaryTableDF,by="ID",all=TRUE)
  my.na <- is.na(summaryTableDFFinal$Freq.y)
  summaryTableDFFinal$Freq.y[my.na] <- summaryTableDFFinal$Freq.x[my.na]
  thisNumTP = summaryTableDFFinal$Freq.y

  #hack for ca testing
  #classIdsList = c(1,2,3,4,5,6)
  #thisNumTP = c(82,99,100,100,100,100)

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

  #Stdev of weights calculation for Wplus. P321 in Bonham-Carter
  VarianceWplus = (1.0 / thisNumTP) + (1.0 / (thisNumAreaUnits-thisNumTP))
  SWplus = sqrt(VarianceWplus)
  SWplus

  #Stdev of weights calculation for Wminus. P321 in Bonham-Carter
  VarianceWminus = (1.0 / (NumTP-thisNumTP)) + (1.0 / (TotalArea-thisNumAreaUnits-NumTP+thisNumTP))
  SWminus = sqrt(VarianceWminus)
  SWminus

  #Contrast vector
  Contrast = Wplus - Wminus

  #Contrast Stdev S(C) i.e. Standard deviation of contrast S(C) = S(var(W+) + var(W-))
  SConstrast = sqrt(VarianceWplus + VarianceWminus)

  #Studentized Contrast C/S(C)
  StudentizedConstrast = Contrast/SConstrast
  StudentizedConstrast

  #generate a table of these values

  weightsTable = data.frame(cbind(thisNumAreaUnits,thisNumTP , Wplus,  SWplus, Wminus, SWminus,Contrast,StudentizedConstrast))
  weightsTable
  plot(weightsTable$Contrast)


  #Finished determining Determine new class values for Raster
  newClassValues = Wplus #This needs some checking based on a cumululative ascending / descending param

  revisedClassesDf = data.frame(cbind(classIdsList,newClassValues))
  names(revisedClassesDf) = c("id","v")

  #subsitute in the weight values into the class values
  weightedClassRaster <- subs(r, revisedClassesDf)

  return(weightedClassRaster)
}
