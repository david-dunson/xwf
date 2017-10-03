#' Default psi list
#' 
#' List with the same local feature functions psi as in the original paper
#' 
#' @return List with 4 different local features psi
#' 
#' @examples
#' default_psi()
#' 
#' @export
default_psi <-
function() list(
  function(x, t) rep(1, length(t)),
  function(x, t) x(t),
  function(x, t) pmax(x(t)-x(t-1), 0),
  function(x, t) -pmin(x(t)-x(t-1), 0)
)
