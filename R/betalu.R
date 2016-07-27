#' The General Support Beta Distribution
#' 
#' Density, distribution function, quantile function and random 
#' generation for the general support beta distribution.
#' 
#' The general support beta distribution with parameters shape1, 
#' shape2, l, and u is the distribution of the random variable 
#' Y = (u-l)*X + l where X ~ Beta(shape1, shape2).  For details 
#' on that distribution, see \code{\link{dbeta}}.
#' 
#' 
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If length(n) > 1, the length is 
#'   taken to be the number required.
#' @param shape1 beta parameter shape1
#' @param shape2 beta parameter shape2
#' @param l lower bound of support
#' @param u upper bound of support
#' @param ncp	non-centrality parameter.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are 
#'   P[X <= x] otherwise, P[X > x].
#' @importFrom stats dbeta pbeta qbeta rbeta
#' @author David Kahle \email{david.kahle@@gmail.com}
#' @name betalu
#' @examples
#' 
#' s <- seq(-1, 1, .01)
#' plot(s, dbetalu(s, 5, 5, -1, 1), type = "l")
#' 
#' f <- function(x) dbetalu(x, 5, 5, -1, 1)
#' x <- -0.5
#' integrate(f, -1, x)
#' (p <- pbetalu(x, 5, 5, -1, 1))
#' qbetalu(p, 5, 5, -1, 1)
#' mean(rbetalu(1e6, 5, 5, -1, 1) <= -.5) # ~= p
#' 
#' 
#' 
NULL




#' @rdname betalu
#' @export 
dbetalu <- function(x, shape1, shape2, l = 0, u = 1, ncp = 0, log = FALSE) {
  log_f <- dbeta((x-l)/(u-l), shape1, shape2, ncp, log = TRUE) - log(u-l)
  if(log) return(log_f)
  exp(log_f)
}


#' @rdname betalu
#' @export 
pbetalu <- function(q, shape1, shape2, l = 0, u = 1, ncp = 0, lower.tail = TRUE) {
  pbeta((q-l)/(u-l), shape1, shape2, ncp, lower.tail = lower.tail)
}



#' @rdname betalu
#' @export 
qbetalu <- function(p, shape1, shape2, l = 0, u = 1, ncp = 0) {
  (u-l)*qbeta(p, shape1, shape2, ncp) + l
}



#' @rdname betalu
#' @export 
rbetalu <- function(n, shape1, shape2, l = 0, u = 1, ncp = 0) {
  (u-l)*rbeta(n, shape1, shape2, ncp) + l
}

