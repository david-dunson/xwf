#' p-value computation for XWFs
#' 
#' Randomization method to compute p-values for an optimized extrema-weighted features generalized additive model fit.
#' 
#' @param GAMobject The GAMobject returned by \code{\link{xwfGridsearch}}
#' @param xx List of function for which to compute the XWFs
#' @param t Matrix containing the times at which the functions xx were measured: Element (i,j) contains the time of the j-th measurement of the i-th function.
#' @param n.i Vector containing the number of measurements for each function. The first n.i[i] elements of the i-th row of t should not be NA.
#' @param psi.list List of predefined local features which are functions of a function (first argument) and a measurement time (second argument)
#' @param F CDF of the values of the functions xx
#' @param n.boot Number for randomizations used to obtain the p-values. The resolution of the p-values is 1/n.boot
#' @param progressbar Boolean specifying whether a progress bar indicating which randomizations have been completed should be displayed.
#' 
#' @return Named vector with p-values
#' 
#' @examples
#' # Data simulation similar to Section 3.2 of the paper
#' 
#' # Sample size
#' n <- 100
#' 
#' # Length of trajectories
#' n.i <- rep(5, n)
#' max.n.i <- max(n.i)
#' 
#' # Times
#' t <- matrix(NA_integer_, nrow = n, ncol = max.n.i)
#' for(i in 1:n) t[i, 1:n.i[i]] <- 1:n.i[i]
#' 
#' 
#' # Sample periods
#' phi <- runif(n = n, min = 1, max = 10)
#' 
#' # Sample offsets
#' m <- 10*runif(n = n)
#' 
#' # Blood pressure measurements
#' x <- t
#' for(i in 1:n) x[i, 1:n.i[i]] <- sin(phi[i] * 2*pi/max.n.i * t[i, 1:n.i[i]]) + m[i]
#' 
#' # Matrix with covariates z
#' q <- 2 # Number of covariates
#' z <- matrix(rnorm(n = n*q), nrow = n, ncol = q)
#' 
#' # Generate outcomes
#' temp <- phi*min(m, 7)
#' temp <- 40*temp
#' prob <- 1/(1+exp( 2*( median(temp)-temp ) ))
#' y <- rbinom(n = n, size = 1, prob = prob)
#' 
#' xx <- list()
#' for(i in 1:n) xx[[i]] <- approxfun(x = t[i,1:n.i[i]], y = x[i,1:n.i[i]], rule = 2)
#' 
#' # Estimate f
#' weights <- matrix(1/n.i, ncol = max.n.i, nrow = n)[!is.na(t)]
#' f <- density(
#' x = t(sapply(X = 1:n, FUN = function(i) c(xx[[i]](t[i,1:n.i[i]]), rep(NA, max.n.i-n.i[i])))),
#' weights = weights/sum(weights),
#' na.rm = T
#' )
#' 
#' # Define CDF of f, F
#' CDF <- c(0)
#' for(i in 2:length(f$x)) CDF[i] <- CDF[i-1]+(f$x[i]-f$x[i-1])*(f$y[i]+f$y[i-1])/2
#' F <- approxfun(x = f$x, y = CDF/max(CDF), yleft = 0, yright = 1)
#' 
#' psi <- list(
#'   function(x, t) abs(x(t)-x(t-1))
#' )
#' 
#' XWFresult <- xwfGridsearch(y = y, xx = xx, t = t, n.i = n.i, psi.list = psi, F = F)
#' 
#' \donttest{XWFpValues(
#' GAMobject = XWFresult$GAMobject,
#' xx = xx,
#' t = t,
#' n.i = n.i,
#' psi.list = psi,
#' F = F,
#' n.boot = 3
#' )}
#' 
#' @importFrom utils txtProgressBar setTxtProgressBar
#' 
#' @export
XWFpValues <- function(GAMobject, xx, t, n.i, F, psi.list = NULL, n.boot = 100, progressbar = TRUE) {
  if(!is.function(F)) stop("marginal CDF 'F' is not specified as a function")
  
  if(is.null(psi.list)) {
    warning("'psi.list' is not specified. Therefore using default.psi(). Make sure psi.list are the same as used to obtain 'GAMobject' for the p-valus to be correct.")
    psi.list <- default_psi()
  }
  
  y <- GAMobject$model[,1]
  
  n.p <- length(GAMobject$model)
  pValMat <- matrix(NA_real_, nrow = n.boot, ncol = n.p)
  
  if(progressbar) pb <- txtProgressBar(max = n.boot, style = 3)
  for(s in 1:n.boot) {
    
    # We simply randomize y to get bootstrap samples
    temp <- xwfGridsearch(y = sample(y), xx = xx, t = t, n.i = n.i, psi.list = psi.list, F = F, progressbar = FALSE)$GAMobject
    temp <- c(summary(temp)$p.pv, summary(temp)$s.pv)
    if(n.p == length(temp)) pValMat[s, ] <- temp # This check on 'temp' is here to ensure 'xwfGridsearch did not accidentially return NA
    
    if(progressbar) setTxtProgressBar(pb, s)
  }
  
  if(progressbar) close(pb)
  
  pValues <- rep(NA_real_, n.p)
  names(pValues) <- c(names(summary(GAMobject)$p.pv), rownames(summary(GAMobject)$s.table))
  
  temp.p <- c(summary(GAMobject)$p.pv, summary(GAMobject)$s.pv)
  for(s in 1:n.p) {
    temp <- temp.p[s] >= pValMat[, s]
    pValues[s] <- (sum(temp, na.rm = TRUE)+1)/sum(!is.na(temp))
  }
  
  return(pValues)
}