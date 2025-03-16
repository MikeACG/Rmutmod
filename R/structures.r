new_MutMatrix <- function(
    modeldt = data.table::data.table(),
    mafdir = character(1L),
    cohort = character(0L),
    k = integer(1L),
    targetdir = character(1L),
    genomedir = character(1L),
    chrs = character(0L),
    fdirs = setNames(character(0L), character(0L)),
    fplabs = setNames(list(), character(0L))
) {

    mutMatrix <- structure(
        list(
            modeldt = modeldt,
            mafdir = mafdir,
            cohort = cohort,
            k = k,
            targetdir = targetdir,
            genomedir = genomedir,
            chrs = chrs,
            fdirs = fdirs,
            fplabs = fplabs
        ),
        class = c("Rmutmod", "MutMatrix")
    )

    return(mutMatrix)

}

validate_MutMatrix <- function(mutMatrix) {

    if (!("data.table" %in% class(mutMatrix$modeldt))) stop("model must be a data.table")
    .mcols <- c("kmer", "ref", "mut", "n", "abundance.adj", "mutRate")

    # check mandatory character columns
    if (!(.mcols[1:3] %in% names(mutMatrix$modeldt))) stop("missing columns in model")
    if (!(is.character(mutMatrix$modeldt$kmer) & is.character(mutMatrix$modeldt$ref) & is.character(mutMatrix$modeldt$mut))) {

        stop("model columns of wrong type")

    }

    # check mandatory number columns
    if (!(.mcols[4:6] %in% names(mutMatrix$modeldt))) stop("missing columns in model")
    if (!(is.integer(mutMatrix$modeldt$n) & is.numeric(mutMatrix$modeldt$abundance.adj) & is.numeric(mutMatrix$modeldt$mutRate))) {

        stop("model columns of wrong type")

    }

    # check optional character feature columns
    .ocols <- setdiff(names(mutMatrix$modeldt), .mcols)
    if (length(.ocols) > 0) {

        if (!(all(grepl("^f[.]", .ocols)))) stop("optional columns in model incorrectly formatted")
        if (!(all(is.character(.ocols)))) stop("model optional columns must be of type character")

    }

}

new_MultiMAFglmmTMB <- function(
    modelPaths = character(6),
    mafdir = character(1),
    k = integer(1),
    targetdir = character(1),
    genomedir = character(1),
    chrs = character(0),
    fdirs = setNames(character(0), character(0)),
    cohort = character(1)
) {

    multiMAFglmmTMB <- structure(
        list(
            modelPaths = modelPaths,
            mafdir = mafdir,
            k = k,
            targetdir = targetdir,
            genomedir = genomedir,
            chrs = chrs,
            fdirs = fdirs,
            cohort = cohort
        ),
        class = c("Rmutmod", "MultiMAFglmmTMB")
    )

    return(multiMAFglmmTMB)

}

new_MonoMAFglmmTMB <- function(
    model = structure(list(), class = "glmmTMB"),
    mafdir = character(1),
    k = integer(1),
    targetdir = character(1),
    genomedir = character(1),
    chrs = character(0),
    fdirs = setNames(character(0), character(0)),
    cohort = character(1)
) {

    monoMAFglmmTMB <- structure(
        list(
            model = model,
            mafdir = mafdir,
            k = k,
            targetdir = targetdir,
            genomedir = genomedir,
            chrs = chrs,
            fdirs = fdirs,
            cohort = cohort
        ),
        class = c("Rmutmod", "MonoMAFglmmTMB")
    )

    return(monoMAFglmmTMB)

}

new_MonoMAFglmmTMBsim <- function(
    fixef = matrix(0.0, nrow = 0, ncol = 0),
    ranef = list(list()),
    sigma = numeric(1),
    cformula = formula(),
    iformulas = list(),
    rformulas = list(),
    flevels = list()
) {

    monoMAFglmmTMBsim <- structure(
        list(
            "fixef" = fixef,
            "ranef" = ranef,
            "sigma" = sigma,
            "cformula" = cformula,
            "iformulas" = iformulas,
            "rformulas" = rformulas,
            "flevels" = flevels
        ),
        class = c("MonoMAFglmmTMBsim")
    )

    return(monoMAFglmmTMBsim)

}

new_MultiMAFglmmTMBsim <- function(
    sims = rep(new_MonoMAFglmmTMBsim(), 6),
    nsims = integer(1)
) {

    multiMAFglmmTMBsim <- structure(
        list(
            "sims" = sims,
            "nsims" = nsims
        ),
        class = c("Rmutsim", "MultiMAFglmmTMBsim")
    )

    return(multiMAFglmmTMBsim)

}


#' @export
modelPathsGet <- function(x) {

    UseMethod("modelPathsGet", x)

}

#' @export
modelPathsGet.MultiMAFglmmTMB <- function(multiMAFglmmTMB) {

    return(multiMAFglmmTMB$modelPaths)

}

#' @export
kGet <- function(x) {

    UseMethod("kGet", x)

}

#' @export
kGet.Rmutmod <- function(rmutmod) {

    return(rmutmod$k)

}

#' @export
targetdirGet <- function(x) {

    UseMethod("targetdirGet", x)

}

#' @export
targetdirGet.Rmutmod <- function(rmutmod) {

    return(rmutmod$targetdir)

}

#' @export
modelGet <- function(x) {

    UseMethod("modelGet", x)

}

#' @export
modelGet.MutMatrix <- function(mutmatrix) {

    return(mutmatrix$modeldt)

}

#' @export
modelGet.MonoMAFglmmTMB <- function(monoMAFglmmTMB) {

    return(monoMAFglmmTMB$model)

}

#' @export
mutpredict <- function (x, newdata, ...) {

   UseMethod("mutpredict", x)

}

#' @export
mutpredict.MutMatrix <- function(mutmatrix, newdata, ...) {

    modeldt <- modelGet(mutmatrix)
    fnames <- names(fdirsGet(mutmatrix))

    newdata[modeldt, "mutRate" := i.mutRate, on = c("kmer", "mut", fnames)]

    return()

}

#' @export
pkmersGet <- function(x) {

    UseMethod("pkmersGet", x)

}

#' @export
pkmersGet.MutMatrix <- function(mutmatrix) {

    modeldt <- modelGet(mutmatrix)
    return(unique(modeldt$kmer))

}

#' @export
genomedirGet <- function(x) {

    UseMethod("genomedirGet", x)

}

#' @export
genomedirGet.Rmutmod <- function(rmutmod) {

    return(rmutmod$genomedir)

}

#' @export
mafdirGet <- function(x) {

    UseMethod("mafdirGet", x)

}

#' @export
mafdirGet.Rmutmod <- function(rmutmod) {

    return(rmutmod$mafdir)

}

#' @export
cohortGet <- function(x) {

    UseMethod("cohortGet", x)

}

#' @export
cohortGet.Rmutmod <- function(rmutmod) {

    return(rmutmod$cohort)

}

#' @export
chrsGet <- function(x) {

    UseMethod("chrsGet", x)

}

#' @export
chrsGet.Rmutmod <- function(rmutmod) {

    return(rmutmod$chrs)

}


#' @export
fdirsGet <- function(x) {

    UseMethod("fdirsGet", x)

}

#' @export
fdirsGet.Rmutmod <- function(rmutmod) {

    return(rmutmod$fdirs)

}

