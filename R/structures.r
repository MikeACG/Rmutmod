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

new_MutGLMs <- function(
    models = list(),
    mafdir = character(1L),
    cohort = character(1L),
    k = integer(1L),
    targetdir = character(1L),
    genomedir = character(1L),
    chrs = character(0L),
    fdirs = setNames(character(0L), character(0L)),
    fplabs = setNames(list(), character(0L)),
    formula = as.formula(NULL),
    warns = setNames(character(length(models)), character(length(models)))
) {

    mutGLMs <- structure(
        list(
            models = models,
            mafdir = mafdir,
            cohort = cohort,
            k = k,
            targetdir = targetdir,
            genomedir = genomedir,
            chrs = chrs,
            fdirs = fdirs,
            fplabs = fplabs,
            formula = formula,
            warns = warns
        ),
        class = c("Rmutmod", "MutGLMs")
    )

    return(mutGLMs)

}

new_MutGLMMTMB <- function(
    model = structure(list(), class = "glmmTMB"),
    mafdir = character(1L),
    cohort = character(1L),
    k = integer(1L),
    targetdir = character(1L),
    genomedir = character(1L),
    chrs = character(0L),
    fdirs = setNames(character(0L), character(0L)),
    fplabs = setNames(list(), character(0L)),
    formula = as.formula(NULL)
) {

    mutGLMMTMB <- structure(
        list(
            model = model,
            mafdir = mafdir,
            cohort = cohort,
            k = k,
            targetdir = targetdir,
            genomedir = genomedir,
            chrs = chrs,
            fdirs = fdirs,
            fplabs = fplabs,
            formula = formula
        ),
        class = c("Rmutmod", "MutGLMMTMB")
    )

    return(mutGLMMTMB)

}

#' @export
kGet <- function(x) {

    UseMethod("kGet")

}

#' @export
kGet.Rmutmod <- function(rmutmod) {

    return(rmutmod$k)

}

#' @export
targetdirGet <- function(x) {

    UseMethod("targetdirGet")

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
modelGet.MutGLMs <- function(mutGLMs) {

    return(mutGLMs$models)

}

#' @export
modelGet.MutGLMMTMB <- function(mutGLMMTMB) {

    return(mutGLMMTMB$model)

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

#' @export
fplabsGet <- function(x) {

    UseMethod("fplabsGet", x)

}

#' @export
fplabsGet.Rmutmod <- function(rmutmod) {

    return(rmutmod$fplabs)

}

#' @export
formulaGet <- function(x) {

    UseMethod("formulaGet", x)

}

#' @export
formulaGet.MutGLMs <- function(mutGLMs) {

    return(mutGLMs$formula)

}

#' @export
formulaGet.MutGLMMTMB <- function(mutGLMMTMB) {

    return(mutGLMMTMB$formula)

}

#' @export
nparamGet <- function(x) {

    UseMethod("nparamGet", x)

}

#' @export
nparamGet.MutMatrix <- function(mutMatrix) {

    modeldt <- modelGet(mutMatrix)
    
    return(nrow(modeldt) - 1)

}

#' @export
nparamGet.MutGLMs <- function(mutGLMs) {

    models <- modelGet(mutGLMs)
    nparam <- sapply(models, function(m) sum(!is.na(coef(m))))

    return(sum(nparam))

}

#' @export
pkmersGet.MutGLMs <- function(mutGLMs) {

    models <- modelGet(mutGLMs)
    pkmers <- gsub("_.*", "", names(models))

    return(pkmers)

}

#' @export
pmutcatGet <- function(x) {

    UseMethod("pmutcatGet", x)

}

#' @export
pmutcatGet.MutMatrix <- function(mutMatrix) {

    modeldt <- modelGet(mutMatrix)
    catdt <- unique(modeldt[, .SD, .SDcols = c("kmer", "mut")], by = c("kmer", "mut"))

    return(catdt)

}

#' @export
pmutcatGet.MutGLMs <- function(mutGLMs) {

    models <- modelGet(mutGLMs)
    catdt <- setNames(
        data.table::fread(text = names(models), sep = "_", header = FALSE),
        c("kmer", "mut")
    )

    return(catdt)

}
