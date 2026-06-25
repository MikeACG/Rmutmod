#' @export
predict.Rmutmod <- function(rmutmod, xdt) {

    model <- rmutmod$model
    xdt[, "p" := useRmutmod(model, .SD)]

    return()

}

#' @export
useRmutmod <- function(x, ...) {

    UseMethod("useRmutmod", x)

}

#' @export
useRmutmod.RmutregList <- function(rmutregList, .xdt) {

    tmpdt <- data.table::copy(.xdt)
    tmpdt[
        ,
        "p" := useRmutreg(rmutregList[[paste0(.BY$kmer, "_", .BY$mut)]], .SD),
        by = c("kmer", "mut")
    ]

    return(tmpdt$p)

}

#' @export
useRmutreg <- function(rmutreg, mcatdt) {

    .formula <- formulaGet(rmutreg)
    beta_hat <- coef(rmutreg)
    sigma_hat <- sigma(rmutreg)

    # parse formula
    formula_vec <- as.character(.formula)
    formulaString <- paste0("~", gsub("[[:space:]]", "", formula_vec[3]))

    X <- model.matrix(as.formula(formulaString), mcatdt)
    mu <- exp(X %*% beta_hat)
    p <- (mu * sigma_hat) / (mu + sigma_hat)

    return(p[, 1])

}

