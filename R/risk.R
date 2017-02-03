##
## Functions for computing second moments,
## partial derivatives, and contributions
##
#' Second centered moments
#'
#' These functions relate to the computation of second centered moments (\code{M2()}),
#' the portfolio variance risk (\code{pm2()} and \code{PortRisk()}),
#' the partial derivatives (\code{dm2()} and \code{PortRiskDeriv()}) and
#' the risk contributions of the assets (\code{cm2()} and \code{PortRiskContrib()}).
#'
#' @param r \code{matrix}, a \eqn{(T x N)} array of returns.
#' @param w \code{numeric}, a \eqn{(N x 1)} vector of portfolio weights.
#' @param percentage \code{logical}, whether risk contributions are expressed as percentages.
#'
#' @return \code{numeric}
#'
#' @name risk
#' @family risk
#'
#' @references Boudt, K. and Peterson, B. and Croux, C. (2008/09), Estimation and decomposition of downside risk for portfolios with non-normal returns, \emph{The Journal of Risk}, \bold{11}(2), Winter 2008/09, 79--103.
#'
#' Jondeau, E. and Rockinger, M. (2006), Optimal portfolio allocation under higher moments, \emph{European Financial Management}, \bold{12}(1), 29--55.
#'
#' @examples
#' data(EuStockMarkets)
#' r <- diff(log(EuStockMarkets)) * 100
#' N <- ncol(r)
#' w <- rep(1 / N, N) ## equal weight allocation
#' M2(r)
#' pm2(r, w)
#' dm2(r, w)
NULL

#' @rdname risk
#' @export
M2 <- function(r){
    L <- nrow(r)
    rc <- apply(r, 2, scale, scale = FALSE) ## centering
    ans <- 1 / (L - 1) * crossprod(rc)
    ans
}
#' @rdname risk
#' @export
pm2 <- function(r, w){
    c(crossprod(w, M2(r)) %*% w)
}
#' @rdname risk
#' @export
dm2 <- function(r, w){
    ans <- 2 * (M2(r) %*% w)
    colnames(ans) <- "DerivMom2"
    rownames(ans) <- colnames(r)
    ans
}
#' @rdname risk
#' @export
cm2 <- function(r, w, percentage = TRUE){
    if (percentage){
        ans <- w * dm2(r, w) / pm2(r, w)
    } else {
        ans <- w * dm2(r, w)
    }
    colnames(ans) <- "ContribRisk"
    ans / 2.0
}
#' @rdname risk
#' @export
PortRisk <- function(r, w){
    pm2(r, w)
}
#' @rdname risk
#' @export
PortRiskDeriv <- function(r, w){
    dm2(r, w)
}
#' @rdname risk
#' @export
PortRiskContrib <- function(r, w, percentage = TRUE){
    if (percentage){
        ans <- w * PortRiskDeriv(r, w) / PortRisk(r, w)
    } else {
        ans <- w * PortRiskDeriv(r, w)
    }
    colnames(ans) <- "PortRiskContrib"
    ans / 2.0
}
