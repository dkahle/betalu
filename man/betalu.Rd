% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/betalu.R
\name{betalu}
\alias{betalu}
\alias{dbetalu}
\alias{pbetalu}
\alias{qbetalu}
\alias{rbetalu}
\title{The General Support Beta Distribution}
\usage{
dbetalu(x, shape1, shape2, l = 0, u = 1, ncp = 0, log = FALSE)

pbetalu(q, shape1, shape2, l = 0, u = 1, ncp = 0, lower.tail = TRUE)

qbetalu(p, shape1, shape2, l = 0, u = 1, ncp = 0)

rbetalu(n, shape1, shape2, l = 0, u = 1, ncp = 0)
}
\arguments{
\item{x, q}{vector of quantiles.}

\item{shape1}{beta parameter shape1}

\item{shape2}{beta parameter shape2}

\item{l}{lower bound of support}

\item{u}{upper bound of support}

\item{ncp}{non-centrality parameter.}

\item{log}{logical; if TRUE, probabilities p are given as log(p).}

\item{lower.tail}{logical; if TRUE (default), probabilities are 
P[X <= x] otherwise, P[X > x].}

\item{p}{vector of probabilities.}

\item{n}{number of observations. If length(n) > 1, the length is 
taken to be the number required.}
}
\description{
Density, distribution function, quantile function and random 
generation for the general support beta distribution.
}
\details{
The general support beta distribution with parameters shape1, 
shape2, l, and u is the distribution of the random variable 
Y = (u-l)*X + l where X ~ Beta(shape1, shape2).  For details 
on that distribution, see \code{\link{dbeta}}.
}
\examples{

s <- seq(-1, 1, .01)
plot(s, dbetalu(s, 5, 5, -1, 1), type = "l")

f <- function(x) dbetalu(x, 5, 5, -1, 1)
x <- -0.5
integrate(f, -1, x)
(p <- pbetalu(x, 5, 5, -1, 1))
qbetalu(p, 5, 5, -1, 1)
mean(rbetalu(1e6, 5, 5, -1, 1) <= -.5) # ~= p



}
\author{
David Kahle \email{david.kahle@gmail.com}
}

