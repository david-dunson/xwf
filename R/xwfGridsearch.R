#' Adaptive grid search
#' 
#' Adaptive grid search to optimize the weighting functions in the extrema-weighted features.
#' 
#' @param y Vector with binary outcomes data
#' @param xx List of functions for which to compute the XWFs
#' @param t Matrix containing the times at which the functions xx were measured: Element (i,j) contains the time of the j-th measurement of the i-th function.
#' @param n.i Vector containing the number of measurements for each function. The first n.i[i] elements of the i-th row of t should not be NA.
#' @param psi.list List of predefined local features which are functions of a function (first argument) and a measurement time (second argument)
#' @param F CDF of the values of the functions xx. Ignored if weighting function w is not the default.
#' @param z Optional matrix with covariates to be included as linear predictors in the generalized additive model
#' @param w Weighting function. The default is the one used in the original paper. See the default for what the roles of its 3 arguments are.
#' @param iter Number of levels in the adaptive grid search. The resolution in b obtained is 2^{-1-iter}.
#' @param rel.shift Optional relative reduction of the integration range to avoid instabilities at the end of the integration ranges. Set to 0 if no such correction is desired.
#' @param progressbar Boolean specifying whether a progress bar indicating what level of the adaptive grid has been completed should be displayed.
#' 
#' @return List containing the final XWFs (wL and wR), the parameters for the optimal weighting functions (b.left and b.right), and the gmcv::gamObject corresponding to the final optimal generalized additive model fit.
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
#' XWFresult <- xwfGridsearch(y = y, xx = xx, t = t, n.i = n.i, psi.list = psi, F = F, z = z)
#' 
#' summary(XWFresult$GAMobject)
#' XWFresult$b.left
#' XWFresult$b.right
#' 
#' 
#' @importFrom utils txtProgressBar setTxtProgressBar
#' 
#' @export
xwfGridsearch <-
function(y, xx, t, n.i, psi.list = default_psi(), F = NULL, z = NULL, iter = 3, w = function(t, i, b, left) ifelse(left, min(1, (1-F(xx[[i]](t)))/(1-b)), min(1, F(xx[[i]](t))/b)), rel.shift = .001, progressbar = TRUE) {
# Function that finds the optimal weighting functions via an adaptive grid search

if(!is.function(psi.list[[1]])) stop("local feature 'psi.list' is not specified as a list of functions")
if(!is.function(xx[[1]])) stop("trajectories 'xx' is not specified as a list of functions")

n <- length(xx)
if(length(n.i) == 1) n.i <- rep(n.i, n)


# Find the begin and end time for each subject
t.min <- apply(X = t, MARGIN = 1, FUN = min, na.rm = T)
t.max <- apply(X = t, MARGIN = 1, FUN = max, na.rm = T)
t.range <- t.max-t.min

p <- length(psi.list) # Number of local features


# Set up matrix to store XWFs
XWFmatL <- matrix(data = NA_real_, nrow = p*2^(iter+1), ncol = n+2)
XWFmatR <- matrix(data = NA_real_, nrow = p*2^(iter+1), ncol = n+2)
countL <- 0
countR <- 0

# Functions that return the full set of, respectively, left and right XWFs, without recomputing that same XWF twice
findXWFl <- function(p, b) {
  temp <- which(XWFmatL[ , 1] == p & XWFmatL[ , 2] == b)
  if(length(temp) == 1) return(XWFmatL[temp, -(1:2)])
  
  countL <<- countL+1
  XWFmatL[countL, 1:2] <<- c(p, b)
  XWFmatL[countL, -(1:2)] <<- xwf(xx = xx, t = t, n.i = n.i, psi = psi.list[[p]], w = function(t, i) w(t, i, b = b, left = TRUE), t.min = t.min, t.max = t.max, t.range = t.range)
  return(XWFmatL[countL, -(1:2)])
}

findXWFr <- function(p, b) {
  temp <- which(XWFmatR[ , 1] == p & XWFmatR[ , 2] == b)
  if(length(temp) == 1) return(XWFmatR[temp, -(1:2)])
  
  countR <<- countR+1
  XWFmatR[countR, 1:2] <<- c(p, b)
  XWFmatR[countR, -(1:2)] <<- xwf(xx = xx, t = t, n.i = n.i, psi = psi.list[[p]], w = function(t, i) w(t, i, b = b, left = FALSE), t.min = t.min, t.max = t.max, t.range = t.range)
  return(XWFmatR[countR, -(1:2)])
}


# Initialize wL and wR
wL <- matrix(nrow = n, ncol = p)
wR <- wL

for(pp in 1:p) {
  wL[, pp] <- findXWFl(p = pp, b = .25)
  wR[, pp] <- findXWFr(p = pp, b = .75)
}

# Compute coresponding GAM UBRE score
modelfit <- xwfGAM(wL = wL, wR = wR, y = y, z = z)
score <- modelfit$gcv.ubre

# Intialize weight function parameters (what we're trying to optimize)
left <- rep(.25, p)
right <- rep(.75, p)

if(progressbar) pb <- txtProgressBar(max = iter, style = 3)
for(s in 1:iter) {
  
  # Width of the grid at the current iteration
  grid.width <- .5/2^s
  
  # Optimize each XWF separately
  for(pp in 1:p) {
    
    # First, the left XWF
    wLl <- wL
    wLl[, pp] <- findXWFl(p = pp, b = left[pp] - grid.width)
    modelfit.l <- tryCatch(
      xwfGAM(wL = wLl, wR = wR, y = y, z = z),
      error = function(e) {
        cat("ERROR :",conditionMessage(e), "\n")
        return(modelfit)
      }
    )
    score.l <- modelfit.l$gcv.ubre
    
    wLr <- wL
    wLr[, pp] <- findXWFl(p = pp, b = left[pp] + grid.width)
    modelfit.r <- tryCatch(
      xwfGAM(wL = wLr, wR = wR, y = y, z = z),
      error = function(e) {
        cat("ERROR :",conditionMessage(e), "\n")
        return(modelfit)
      }
    )
    score.r <- modelfit.l$gcv.ubre
    
    temp <- which.max(c(score.l, score.r, score))
    
    if(temp == 1) {
      left[pp] <- left[pp] - grid.width/2
      wL <- wLl
      score <- score.l
      modelfit <- modelfit.l
    } else if(temp == 2) {
      left[pp] <- left[pp] + grid.width/2
      wL <- wLr
      score <- score.r
      modelfit <- modelfit.r
    }
    
    
    # Second, the right XWF
    wRl <- wR
    wRl[, pp] <- findXWFr(p = pp, b = right[pp] - grid.width)
    modelfit.l <- tryCatch(
      xwfGAM(wL = wL, wR = wRl, y = y, z = z),
      error = function(e) {
        cat("ERROR :",conditionMessage(e), "\n")
        return(modelfit)
      }
    )
    score.l <- modelfit.l$gcv.ubre
    
    wRr <- wR
    wRr[, pp] <- findXWFr(p = pp, b = right[pp] + grid.width)
    modelfit.r <- tryCatch(
      xwfGAM(wL = wL, wR = wRr, y = y, z = z),
      error = function(e) {
        cat("ERROR :",conditionMessage(e), "\n")
        return(modelfit)
      }
    )
    score.r <- modelfit.r$gcv.ubre
    
    temp <- which.max(c(score.l, score.r, score))
    
    if(temp == 1) {
      right[pp] <- right[pp] - grid.width/2
      wR <- wRl
      score <- score.l
      modelfit <- modelfit.l
    } else if(temp == 2) {
      right[pp] <- right[pp] + grid.width/2
      wR <- wRr
      score <- score.r
      modelfit <- modelfit.r
    }
  }
  
  if(progressbar) setTxtProgressBar(pb, s)
}
if(progressbar) close(pb)

return(list(wL = wL, wR = wR, b.left = left, b.right = right, GAMobject = modelfit))
}