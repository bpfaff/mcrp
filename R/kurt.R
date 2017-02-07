##
## Functions for computing fourth moments,
## partial derivatives, and contributions
##
#' Fourth centered moments
#'
#' These functions relate to the computation of the kurtosis (\code{M4()}),
#' the portfolio variance risk (\code{pm4()} and \code{PortKurt()}),
#' the partial derivatives (\code{dm4()} and \code{PortKurtDeriv()}) and
#' the risk contributions of the assets (\code{cm4()} and \code{PortKurtContrib()}).
#'
#' @param r \code{matrix}, a \eqn{(T x N)} array of returns.
#' @param w \code{numeric}, a \eqn{(N x 1)} vector of portfolio weights.
#' @param percentage \code{logical}, whether risk contributions are expressed as percentages.
#'
#' @return \code{numeric}
#'
#' @name kurt
#' @family kurt
#'
#' @references Boudt, K. and Peterson, B. and Croux, C. (2008/09),
#' Estimation and decomposition of downside risk for portfolios with non-normal returns,
#' \emph{The Journal of Risk}, \bold{11}(2), Winter 2008/09, 79--103.
#'
#' Jondeau, E. and Rockinger, M. (2006), Optimal portfolio allocation under higher moments,
#' \emph{European Financial Management}, \bold{12}(1), 29--55.
#'
#' @examples
#' data(MultiAsset)
#' MA <- as.timeSeries(MultiAsset[, 1:4])
#' r <- na.omit(diff(log(MA)) * 100)
#' N <- ncol(r)
#' w <- rep(1 / N, N) ## equal weight allocation
#' M4(r)
#' pm4(r, w)
#' dm4(r, w)
NULL

#' @rdname kurt
#' @export
M4 <- function(r){
    N <- ncol(r)
    L <- nrow(r)
    rc <- apply(r, 2, scale, scale = FALSE) ## centering
    ans <- matrix(0, nrow = N, ncol = N ^ 3)
    for (i in 1:L) {
        term1 <- tcrossprod(rc[i, ])
        term2 <- kronecker(term1, t(rc[i, ]))
        ans <- ans + kronecker(term2, t(rc[i, ]))
    }
    ans <- ans / L
    ans
}
#' @rdname kurt
#' @export
pm4 <- function(r, w){
    c(crossprod(w, M4(r)) %*% kronecker(kronecker(w, w), w))
}
#' @rdname kurt
#' @export
dm4 <- function(r, w){
    ans <- 4 * (M4(r) %*% kronecker(kronecker(w, w), w))
    colnames(ans) <- "DerivMom4"
    rownames(ans) <- colnames(r)
    ans
}
#' @rdname kurt
#' @export
cm4 <- function(r, w, percentage = TRUE){
    if (percentage) {
        ans <- w * dm4(r, w) / pm4(r, w)
    } else {
        ans <- w * dm4(r, w)
    }
    colnames(ans) <- "ContribKurt"
    ans / 4.0
}
#' @rdname kurt
#' @export
PortKurt <- function(r, w){
    nomin <- pm4(r, w)
    denom <- pm2(r, w) ^ 2
    ans <- nomin / denom
    ans
}
#' @rdname kurt
#' @export
PortKurtDeriv <- function(r, w){
    term1 <- pm2(r, w) * dm4(r, w)
    term2 <- pm4(r, w) * dm2(r, w)
    term3 <- 2 * pm2(r, w) ^ 3
    ans <- (term1 - term2) / term3
    colnames(ans) <- "PortKurtDeriv"
    rownames(ans) <- colnames(r)
    ans
}
#' @rdname kurt
#' @export
PortKurtContrib <- function(r, w, percentage = TRUE){
    if (percentage) {
        ans <- w * PortKurtDeriv(r, w) / PortKurt(r, w)
    } else {
        ans <- w * PortKurtDeriv(r, w)
    }
    colnames(ans) <- "PortKurtContrib"
    ans
}
