#' Calc Weights Function
#'
#' This function calculates WoE weights of a raster. Currently doesn't handle missing data.
#' @param evidenceRaster Raster of of evidence
#' @param predictTPts Training points
#' @param isCumulative boolean of whether to accumulate area values (e.g for ascending variables) or not (e.g. categorical vars)
#' @keywords WoE, Weights of Evidence
#' @export
#' @examples
#' calcWeights()
calcWeights = function (evidenceRaster,predictTPts,isCumulative) {
  r = raster(evidenceRaster)
  #total number of TP
  NumOfTP = length(predictTPts)

  #Compute the area for each class
  freqOfRDf = data.frame(freq(r))
  classIdsList = freqOfRDf$value
  thisNumAreaUnits = freqOfRDf$count

  #Determine if doing cumulative ascending or categorical
  if (isCumulative==TRUE) {
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
  if (isCumulative==TRUE) {
    summaryTableDF$Freq = cumsum(summaryTableDF$Freq)
  }


  df = data.frame(classIdsList, rep(NumOfTP, length(classIdsList)))
  colnames(df) = c("ID","Freq")
  colnames(summaryTableDF) = c("ID","Freq")
  summaryTableDFFinal = merge(df,summaryTableDF,by="ID",all=TRUE)
  my.na <- is.na(summaryTableDFFinal$Freq.y)
  summaryTableDFFinal$Freq.y[my.na] <- summaryTableDFFinal$Freq.x[my.na]
  thisNumOfTP = summaryTableDFFinal$Freq.y


  #Trap from ArcSDM to catch situations where no TPs in the class
  logical =  thisNumOfTP == NumOfTP
  thisNumOfTP[logical] = thisNumOfTP[logical] - 0.01

  TotalAreaVector = rep(TotalArea,length(thisNumOfTP))

  #P(D)
  PD = NumOfTP/TotalAreaVector

  #P(!D)
  PNotD = (TotalAreaVector - NumOfTP)/TotalAreaVector

  #Conditional probabilities for W+ and W-
  #####################

  #Calculate the weight+
  #####################
    #P(B|D)
  PB_D = thisNumOfTP/NumOfTP

  #P(B|!D)
  PB_NotD = (thisNumAreaUnits-thisNumOfTP) / (TotalAreaVector-NumOfTP)
  PB_NotD


  #Wplus = LN(P(B|D) / P(B|!D))
  Wplus = log(PB_D/PB_NotD)
  Wplus


  #Calculate the weight-
  ####################
  #P(!B|D)
  PNotB_D = (NumOfTP-thisNumOfTP)/NumOfTP
  PNotB_D

  #P(!B|!D)
  PNotB_NotD = (TotalAreaVector-thisNumAreaUnits-NumOfTP+thisNumOfTP)/(TotalAreaVector-NumOfTP)
  PNotB_NotD

  #Wminus = LN(P(!B|D) / P(!B|!D))
  Wminus = log(PNotB_D/PNotB_NotD)
  Wminus



  #Stdev of weights calculation for Wplus. P321 in Bonham-Carter
  VarianceWplus = (1.0 / thisNumOfTP) + (1.0 / (thisNumAreaUnits-thisNumOfTP))
  SWplus = sqrt(VarianceWplus)
  SWplus

  #Stdev of weights calculation for Wminus. P321 in Bonham-Carter
  VarianceWminus = (1.0 / (NumOfTP-thisNumOfTP)) + (1.0 / (TotalArea-thisNumAreaUnits-NumOfTP+thisNumOfTP))
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
  weightsTable = data.frame(cbind(thisNumAreaUnits,thisNumOfTP , Wplus,  SWplus, Wminus, SWminus,Contrast,StudentizedConstrast))
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
