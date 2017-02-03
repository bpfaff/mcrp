##
## Functions for computing third moments,
## partial derivatives, and contributions
##
#' Third centered moments
#'
#' These functions relate to the computation of (co-)skewness (\code{M3()}),
#' the portfolio skewness (\code{pm3()} and \code{PortSkew()}), the partial
#' derivatives (\code{dm3()} and \code{PortSkewDeriv()})
#' and the contributions of the assets to skewness (\code{cm3()} and
#' \code{PortSkewContrib()}).
#'
#' @param r \code{matrix}, a \eqn{(T x N)} array of returns.
#' @param w \code{numeric}, a \eqn{(N x 1)} vector of portfolio weights.
#' @param percentage \code{logical}, whether risk contributions are expressed as percentages.
#'
#' @return \code{numeric}
#'
#' @name skew
#' @family skew
#'
#' @references Boudt, K. and Peterson, B. and Croux, C. (2008/09), Estimation and decomposition of
#' downside risk for portfolios with non-normal returns, \emph{The Journal of Risk}, \bold{11}(2),
#' Winter 2008/09, 79--103.
#'
#' Jondeau, E. and Rockinger, M. (2006), Optimal portfolio allocation under higher moments,
#' \emph{European Financial Management}, \bold{12}(1), 29--55.
#'
#' @examples
#' data(EuStockMarkets)
#' r <- diff(log(EuStockMarkets)) * 100
#' N <- ncol(r)
#' w <- rep(1 / N, N) ## equal weight allocation
#' M3(r)
#' pm3(r, w)
#' dm3(r, w)
NULL

#' @rdname skew
#' @export
M3 <- function(r){
    N <- ncol(r)
    L <- nrow(r)
    rc <- apply(r, 2, scale, scale = FALSE) ## centering
    ans <- matrix(0, nrow = N, ncol = N ^ 2)
    for (i in 1:L) {
        ans <- ans + kronecker(tcrossprod(rc[i, ]), t(rc[i, ]))
    }
    ans <- ans / L
    ans
}
#' @rdname skew
#' @export
pm3 <- function(r, w){
    ans <- c(crossprod(w, M3(r) %*% kronecker(w, w)))
    ans
}
#' @rdname skew
#' @export
dm3 <- function(r, w){
    ans <- 3 * (M3(r) %*% kronecker(w, w))
    colnames(ans) <- "DerivMom3"
    rownames(ans) <- colnames(r)
    ans
}
#' @rdname skew
#' @export
cm3 <- function(r, w, percentage = TRUE){
    if (percentage) {
        ans <- w * dm3(r, w) / pm3(r, w)
    } else {
        ans <- w * dm3(r, w)
    }
    colnames(ans) <- "ContribSkew"
    ans / 3.0
}
#' @rdname skew
#' @export
PortSkew <- function(r, w){
    nomin <- pm3(r, w)
    denom <- pm2(r, w) ^ (3 / 2)
    nomin / denom
}
#' @rdname skew
#' @export
PortSkewDeriv <- function(r, w){
    term1 <- pm2(r, w) ^ (3 / 2) * dm3(r, w)
    term2 <- pm3(r, w) * pm2(r, w) ^ 0.5 * dm2(r, w)
    term3 <- pm2(r, w) ^ 3
    ans <- (term1 - term2) / term3
    colnames(ans) <- "PortSkewDeriv"
    rownames(ans) <- colnames(r)
    ans
}
#' @rdname skew
#' @export
PortSkewContrib <- function(r, w, percentage = TRUE){
    if (percentage) {
        ans <- w * PortSkewDeriv(r, w) / PortSkew(r, w)
    } else {
        ans <- w * PortSkewDeriv(r, w)
    }
    colnames(ans) <- "PortSkewContrib"
    ans
}
