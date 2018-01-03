#' Evaluate the GAM
#' 
#' Evaluate the generalized additive model for a set of computed extrema-weighted features
#' 
#' @param wL Matrix with left extrema-weighted features
#' @param wR Matrix with right extrema-weighted features
#' @param y Binary vector with outcomes
#' @param z Optional matrix z with extra, linear predictors
#' 
#' @examples
#' xwf:::xwfGAM(wL = rep(1:45, 10), wR = rep(1:90, 5), y = c(rep(0:1, 225)))
#' 
#' @importFrom mgcv gam
xwfGAM <- function(wL, wR, y, z = NULL) {
  p <- ncol(wL)
  if(is.null(p)) p <- 1
  if(!is.null(z)) {
    q <- ncol(z)
    if(is.null(q)) q <- 1
  }
  
  eval(parse(text = paste0(
    ifelse(is.null(z),
      "mgcv::gam(y ~ 1",
      ifelse(q == 1,
        "mgcv::gam(y ~ s(z)",
        paste0(c("mgcv::gam(y ~ s(z[,1])", sapply(X = 2:q, FUN = function(j) paste0("+s(z[,", j, "])"))), collapse = '')
      )
    ),
    ifelse(p == 1,
      "+s(wL)+s(wR)",
      paste0(sapply(X = 1:p, FUN = function(j) paste0("+s(wL[,", j, "])+s(wR[,", j, "])")), collapse = '')
    ),
    ", family = binomial(link = 'logit'))"
  )))
}
