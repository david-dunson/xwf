#' Compute XWFs
#' 
#' Compute extrema-weighted features based on functions, predefined local features, and weighting functions
#' 
#' @param xx List of function for which to compute the XWFs
#' @param t Matrix containing the times at which the functions xx were measured: Element (i,j) contains the time of the j-th measurement of the i-th function.
#' @param n.i Vector containing the number of measurements for each function. The first n.i[i] elements of the i-th row of t should not be NA.
#' @param psi Predefined local feature which is a function of a function (first argument) and a measurement time (second argument)
#' @param w Weighting function. The default is the one used in the original paper.
#' @param b Parameter of the weighting function. See original paper for details. Ignored if weighting function w is not the default.
#' @param F CDF of the values of the functions xx. Ignored if weighting function w is not the default.
#' @param t.min Vector with time of first measurement for each function. Computed from t if omitted but providing it saves computational cost.
#' @param t.max Analogous to t.min but now the time of the last measurement.
#' @param t.range Vector with differences between t.max and t.min. Can be supplied to avoid recomputation.
#' @param rel.shift Optional relative reduction of the integration range to avoid instabilities at the end of the integration ranges. Set to 0 if no such correction is desired.
#' @param left Boolean specifying whether the left (TRUE) or right (FALSE) extrema-weighted features should be computed: Left and right refer to the weighting function. Ignored if weighting function w is not the default.
#' 
#' @return Vector containing the extrema-weighted features obtained by numerical integration for each of the functions.
#' 
#' @examples
#' xwf(
#' xx = list(function(t) t),
#' t = (1:10)/10,
#' n.i = 10,
#' psi = function(x, t) x(t),
#' b = .2,
#' F = function(x) x
#' )
#' 
#' @importFrom stats integrate
#' 
#' @export
xwf <-
function(xx, t, n.i, psi, w = function(t, i) ifelse(left, min(1, (1-F(xx[[i]](t)))/(1-b)), min(1, F(xx[[i]](t))/b)), b = .5, F = NULL, t.min = NULL, t.max = NULL, t.range = NULL, rel.shift = .001, left = TRUE) {
  # This R script contains a function that computes Extrema Weighted Features (XWF)
  
  if(!is.null(F) &!is.function(F)) stop("marginal CDF 'F' is not specified as a function")
  if(!is.function(psi)) stop("local feature 'psi' is not specified as a function")
  if(!is.function(xx[[1]])) stop("trajectories 'xx' is not specified as a list of functions")
  
  n <- length(xx)
  if(length(n.i) == 1) n.i <- rep(n.i, n)
  
  # Compute t.min, t.max, and t.range if they are given.
  if(is.null(dim(t))) {
    if(is.null(t.min)) t.min <- min(t, na.rm = T)
    if(is.null(t.max)) t.max <- max(t, na.rm = T)
  } else {
    if(is.null(t.min)) t.min <- apply(X = t, MARGIN = 1, FUN = min, na.rm = T)
    if(is.null(t.max)) t.max <- apply(X = t, MARGIN = 1, FUN = max, na.rm = T)
  }
  if(is.null(t.range)) t.range <- t.max-t.min
  
  # Define the integral for the XWFs
  shift <- rel.shift*t.range
  # The lower and upper bound are shifted slightly inwards to prevent problems due to psi_j(x, t) not being properly defines, for instance at the boundaries.
  
  # Return the vector with XWFs
  return(Vectorize(
    FUN = function(i) integrate(f = function(t) w(t, i)*psi(xx[[i]], t), lower = t.min[i]+shift[i], upper = t.max[i]-shift[i], subdivisions = 100, rel.tol = .05, abs.tol = .05, stop.on.error = FALSE)$value/t.range[i]
  )(1:n))
}