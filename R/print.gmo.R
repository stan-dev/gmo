#' Print objects of class \code{gmo}.
#'
#' @param x object of class \code{gmo}.
#' @param \dots further arguments passed to or from other methods.
#'
#' @export
print.gmo <- function(x, ...) {
  print(x$par, ...)
}
