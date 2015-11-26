#' Otsu binarization
#'
#' Unfinished and not very effective Otsu binarization script converted from http://www.biomecardio.com/matlab/otsu.m
#' @param grayValueRaster for thresholding
#' @keywords thresholding automatically
#' @export
#' @examples
#' otsuMethod()

otsuMethod <-
  function(grayValueRaster){
    inputRaster = readGDAL("mndwiRasterDOSReflectanceNDVI.tif")
    ir = raster(inputRaster)
    I = ir@data@values

    nbins = 256 #Number of bins

    #if (min(I(:)) == max(I(:)) ), disp('otsu error: no intensity variability'); thresh =min(I(:)); return; end;

    intercept = min(I,na.rm = TRUE); #we will translate min-val to be zero
    slope = (nbins-1)/ (max(I,na.rm = TRUE)-intercept); #we will scale images to range 0..(nbins-1)
    # Convert to 256 levels
    I = round((I - intercept) * slope);

    #Probability distribution
    his = hist(I,256,plot=FALSE,right = TRUE)
    histo= his$counts

    P = histo/sum(histo)
    #histo = c(99, his$counts)


    # Zeroth- and first-order cumulative moments
    w = cumsum(P);
    vec = (1:nbins)
    mu = cumsum(vec*P);

    end <- function(x) { return( x[length(x)] ) }

    sigma2B =(end(mu)*w[2:length(w)-1]-mu[2:length(mu)-1])^2/w[2:length(w)-1]/(1-w[2:length(w)-1])


    sigma2BNEW =(mu[length(mu)]*w[2:length(w)-1]-mu[2:length(mu)-1])^2/w[2:length(w)-1]/(1-w[2:length(w)-1])

    sigma2BNEW =(mu[length(mu)]*w[2:length(w)]-mu[2:length(mu)-1])^2/w[2:length(w)-1]/(1-w[2:length(w)-1])


    sigma2B =(mu(end)*w(2:end-1)-mu(2:end-1)).^2./w(2:end-1)./(1-w(2:end-1));


    ^2/w(2:end-1)/(1-w(2:end-1));


    [maxsig,k] = max(sigma2B);
  }
