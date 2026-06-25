#' @export
coef.Rmutreg <- function(rmutreg) {

    return(rmutreg$beta_hat)
    
}

#' @export
sigma.Rmutreg <- function(rmutreg) {

    return(rmutreg$sigma_hat)
    
}

#' @export
formulaGet <- function(x) {

    UseMethod("formulaGet", x)

}

#' @export
formulaGet.Rmutreg <- function(rmutreg) {

    return(rmutreg$formula)

}
