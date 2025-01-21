# # these are a collection of functions that help fit a glm to all
# # possible mutations in the target, currently they are unused
# # because the model is not very scalable once interactions and
# # more than a few features are included (if the target is for instance the whole exome)

# formatFeatures <- function(xdt, pkmers, nflank, fplabs) {

#     # format mutation category
#     xdt[, "mutcat" := stringi::stri_join(kmer, mut, sep = ">")]

#     # specify levels of variables
#     catdt <- expandMuts(data.table::data.table(kmer = pkmers), nflank)
#     catdt[, "mutcat" := stringi::stri_join(kmer, mut, sep = ">")]
#     fpl <- append(list("mutcat" = catdt$mutcat), fplabs)
#     for (jj in 1:length(fpl)) {

#         cname <- names(fpl)[jj]
#         if (length(fpl[[jj]]) > 0) xdt[, (cname) := factor(get(cname), fpl[[jj]])]

#     }

#     return()

# }

vector2sdt <- function(x, jcum) {

    UseMethod("vector2sdt", x)

}

vector2sdt.factor <- function(x, jcum) {

    # initialize data table
    .sdt <- data.table::data.table("i" = 1:length(x))

    # add dummy column where a 1 should be placed for each observation
    .sdt[, "j" := match(x, levels(x))]

    # remove observations of last dummy column and offset column index
    .sdt <- .sdt[j != length(levels(x))]
    .sdt <- .sdt[, "j" := j + jcum]

    # remaining indices should have a 1 as their values
    .sdt[, "v" := 1L]

    return(.sdt)

}

vector2sdt.default <- function(x, jcum) {

    # initialize data table
    .sdt <- data.table::data.table("i" = 1:length(x))

    # add column where values will go offseted by previous columns
    # also add the values themselves
    .sdt[, ':=' ("j" = jcum + 1L, "v" = x)]

    # remove observations with value of 0
    .sdt <- .sdt[v != 0]

    return(.sdt)

}

interaction2sdt <- function(x, jcum, mutcats) {

    UseMethod("interaction2sdt", x)

}

interaction2sdt.factor <- function(x, jcum, mutcats) {

    # initialize data table
    .sdt <- data.table::data.table("i" = 1:length(x))

    # add in what level each observation is for mutcats
    .sdt[, "jmutcat" := match(mutcats, levels(mutcats))]

    # add in what level each observation is for the target feature
    .sdt[, "jx" := match(x, levels(x))]

    # remove observations from last level of any of mutcats or x
    .sdt <- .sdt[jmutcat != length(levels(mutcats)) & jx != length(levels(x))]

    # jx - 1: "sets" of ncolsMutcats "behind" the set where the 1 should go
    # the position inside the set is given by jmutcat
    ncolsMutcats <- length(levels(mutcats)) - 1
    .sdt[, "j" := jmutcat + (ncolsMutcats * (jx - 1L))]

    # remove temporary columns and offset all column indices by the cumulative of columns
    .sdt[, ':=' ("jmutcat" = NULL, "jx" = NULL)]
    .sdt <- .sdt[, "j" := j + jcum]

    # remaining indices should have a 1 as their values
    .sdt[, "v" := 1L]

    return(.sdt)

}

interaction2sdt.default <- function(x, jcum, mutcats) {

    # initialize data table
    .sdt <- data.table::data.table("i" = 1:length(x))

    # add in what level each observation is for mutcats
    .sdt[, "j" := match(mutcats, levels(mutcats))]

    # add value of each observation for the target feature
    .sdt[, "v" := x]

    # remove observation if in last level of mutcats or a value of 0 for x
    .sdt <- .sdt[j != length(levels(mutcats)) & v != 0]

    # offset all column indices by the cumulative of columns
    .sdt <- .sdt[, "j" := j + jcum]

    return(.sdt)

}

dt2sdt <- function(xdt, mutcatInteractions = TRUE) {

    # calculate the number of columns that each feature will create in the model matrix
    nf <- sapply(xdt, function(x) if (is.factor(x)) length(levels(x)) - 1 else 1)
    nfc <- head(cumsum(c(0, nf)), -1)

    # get sparse representation
    sdt <- mapply(vector2sdt, xdt, nfc, SIMPLIFY = FALSE)
    sdt <- data.table::rbindlist(sdt)

    # deal with interactions if neccessary
    isdt <- data.table::data.table()
    nfi <- 0
    if (mutcatInteractions) {

        mutcats <- xdt$mutcat
        xdt[, "mutcat" := NULL]
        nfi <- nf[names(xdt)] * (length(levels(mutcats)) - 1)
        nfic <- head(cumsum(c(0, nfi)), -1) + sum(nf)

        isdt <- mapply(interaction2sdt, xdt, nfic, MoreArgs = list(mutcats = mutcats), SIMPLIFY = FALSE)
        isdt <- data.table::rbindlist(isdt)

    }

    return(rbind(sdt, isdt))

}

sdt2disk <- function(sdt, nr, outpath) {

    h <- c(
        "%%MatrixMarket matrix coordinate real general",
        paste(nr, 1343, nrow(sdt), collapse = " ")
    )
    writeLines(h, outpath)
    data.table::fwrite(sdt, outpath, append = TRUE, sep = " ", col.names = FALSE)

}

# trainMutGLM <- function(mafdir, cohort, k, targetdir, genomedir, chrs, fdirs, fplabs, .formula) {

#     pkmers <- makePkmers(k)
#     nflank <- (k - 1) / 2
#     mafdb <- arrow::open_dataset(mafdir)
#     targetdb <- arrow::open_dataset(targetdir)
#     genomePaths <- paste0(genomedir, chrs, ".fasta")
    
#     y <- list()
#     X <- list()
#     for (ii in 1:length(chrs)) {

#         cat(ii, "/", length(chrs), "...\n", sep = "")

#         xdt <- Rmutmod:::chrom2table(chrs[ii], mafdb, cohort, targetdb, genomePaths[ii], nflank, pkmers, fdirs, fplabs)
#         y[[ii]] <- xdt$nmut
#         xdt[, "nmut" := NULL]

#         sdt <- dt2sdt(xdt)
#         nr <- nrow(xdt)
#         rm(xdt)

#         X[[ii]] <- Matrix::sparseMatrix(
#             sdt$i,
#             sdt$j,
#             x = sdt$v,
#             dims = c(nr, 479)
#         )
#         rm(sdt); gc()

#     }
#     y <- do.call(c, y)
#     X <- do.call(rbind, X)
    
#     # fit model (this is very memory intensive for site-level models)
#     glmMut <- glmnet::glmnet(
#         X,
#         y,
#         "poisson",
#         lambda = 0,
#         standardize = FALSE,
#         intercept = TRUE,
#         trace.it = 1
#     )

#     mutGLM <- new_MutGLM(
#         model = glmMut,
#         mafdir = mafdir,
#         cohort = cohort,
#         k = k,
#         targetdir = targetdir,
#         genomedir = genomedir,
#         chrs = chrs,
#         fdirs = fdirs,
#         fplabs = fplabs,
#         formula = .formula
#     )

#     return(mutGLM)

# }

# mutpredict.MutGLM <- function(mutGLM, newdata, ...) {

#     model <- modelGet(mutGLM)
#     k <- kGet(mutGLM)
#     fplabs <- fplabsGet(mutGLM)
#     .formula <- as.formula(paste0("~", paste(labels(terms(formulaGet(mutGLM))), collapse = "+")))

#     pkmers <- makePkmers(k)
#     nflank <- floor(k / 2)
#     formatFeatures(newdata, pkmers, nflank, fplabs)
#     X <- MatrixModels::model.Matrix(.formula, newdata, sparse = TRUE)
#     newdata[, "mutRate" := predict(model, X, type = "response")]

#     return()

# }

# trainMutGLM <- function(mafdir, cohort, k, targetdir, genomedir, chrs, fdirs, fplabs, .formula) {

#     pkmers <- Rmutmod:::makePkmers(k)
#     nflank <- (k - 1) / 2
#     mafdb <- arrow::open_dataset(mafdir)
#     targetdb <- arrow::open_dataset(targetdir)
    
#     genomePaths <- paste0(genomedir, chrs, ".fasta")
#     catdt <- Rmutmod:::expandMuts(data.table::data.table("kmer" = pkmers), nflank)
#     catdt[, "mutcat" := paste(kmer, mut, sep = ">")]
#     .fplabs <- append(fplabs,  list("mutcat" = catdt$mutcat))

#     y <- list()
#     X <- list()
#     for (ii in 1:length(chrs)) {

#         cat(ii, "/", length(chrs), "...\n", sep = "")

#         xdt <- Rmutmod:::chrom2table(chrs[ii], mafdb, cohort, targetdb, genomePaths[ii], nflank, pkmers, fdirs, fplabs)
#         xdt <- xdt[complete.cases(xdt)]
#         y[[ii]] <- xdt$nmut
#         xdt[, ':=' ("start" = NULL, "rangeid" = NULL, "ref" = NULL, "nmut" = NULL)]

#         xdt[, "mutcat" := stringi::stri_join(kmer, mut, sep = ">")]
#         xdt[, ':=' ("kmer" = NULL, "mut" = NULL)]
#         Rmutmod:::formatFeatures(xdt, .fplabs)
#         sdt <- dt2sdt(xdt)

#         X[[ii]] <- spam::spam(
#             sdt,
#             nrow(xdt),
#             1343
#         )
#         rm(xdt, sdt); gc()

#     }

#     .X <- X[[1]]
#     X[[1]] <- NULL
#     ii <- 2
#     while(length(X) > 0) {

#         print(ii)
#         .X <- rbind(.X, X[[1]])
#         X[[1]] <- NULL
#         ii <- ii + 1

#     }

#     X <- do.call(rbind, X)
#     y <- unlist(y, recursive = FALSE, use.names = FALSE)
#     glmMut <- glmnet::glmnet(
#         X[[22]],
#         y[[22]],
#         "poisson",
#         lambda = 0,
#         standardize = FALSE,
#         intercept = TRUE,
#         trace.it = 1
#     )

#     mutGLMs <- new_MutGLMs(
#         models = models,
#         mafdir = mafdir,
#         cohort = cohort,
#         k = k,
#         targetdir = targetdir,
#         genomedir = genomedir,
#         chrs = chrs,
#         fdirs = fdirs,
#         fplabs = fplabs,
#         formula = .formula,
#         warns = warns
#     )

#     return(mutGLMs)

# }

# #' @export
# trainMutCGLM <- function(mafdir, cohort, k, targetdir, genomedir, chrs, fdirs, fplabs, .formula) {

#     pkmers <- Rmutmod:::makePkmers(k)
#     nflank <- (k - 1) / 2
#     mafdb <- arrow::open_dataset(mafdir)
#     targetdb <- arrow::open_dataset(targetdir)
#     genomePaths <- paste0(genomedir, chrs, ".fasta")

#     ccols <- c("kmer", "mut", names(fplabs))
#     cxdt <- list()
#     for (ii in 1:length(chrs)) {

#         cat(ii, "/", length(chrs), "...\n", sep = "")

#         xdt <- Rmutmod:::chrom2table(chrs[ii], mafdb, cohort, targetdb, genomePaths[ii], nflank, pkmers, fdirs, fplabs)
#         cxdt[[ii]] <- xdt[complete.cases(xdt), list("nmut" = sum(nmut), "nchance" = .N), by = ccols]
#         rm(xdt)

#     }
#     cxdt <- data.table::rbindlist(cxdt)[, list("nmut" = sum(nmut), "nchance" = sum(nchance)), by = ccols]

#     cxdt[, "mutcat" := stringi::stri_join(kmer, mut, sep = ">")]
#     cxdt[, ':=' ("kmer" = NULL, "mut" = NULL)]
#     catdt <- Rmutmod:::expandMuts(data.table::data.table("kmer" = pkmers), nflank)
#     catdt[, "mutcat" := paste(kmer, mut, sep = ">")]
#     Rmutmod:::formatFeatures(cxdt, append(fplabs,  list("mutcat" = catdt$mutcat)))

#     cxdt[, "logNchance" := log(nchance)]
#     cxdt[, "nchance" := NULL]
#     formula_string <- paste0("nmut ~ (mutcat * (", paste(names(fplabs), collapse = " + "), ")) + offset(logNchance)")
#     model <- MASS::glm.nb(
#         as.formula(formula_string),
#         cxdt,
#         control = glm.control(trace = TRUE)
#     )

#     mutGLMs <- new_MutGLMs(
#         models = models,
#         mafdir = mafdir,
#         cohort = cohort,
#         k = k,
#         targetdir = targetdir,
#         genomedir = genomedir,
#         chrs = chrs,
#         fdirs = fdirs,
#         fplabs = fplabs,
#         formula = .formula,
#         warns = warns
#     )

#     return(mutGLMs)

# }


# #' @export
# trainMutGLM <- function(mafdir, cohort, k, targetdir, genomedir, chrs, fdirs, fplabs, .formula) {

#     pkmers <- Rmutmod:::makePkmers(k)
#     nflank <- (k - 1) / 2
#     mafdb <- arrow::open_dataset(mafdir)
#     targetdb <- arrow::open_dataset(targetdir)
    
#     genomePaths <- paste0(genomedir, chrs, ".fasta")
#     catdt <- Rmutmod:::expandMuts(data.table::data.table("kmer" = pkmers), nflank)
#     catdt[, "mutcat" := paste(kmer, mut, sep = ">")]
#     .fplabs <- append(fplabs,  list("mutcat" = catdt$mutcat))

#     d <- paste0("./", basename(tempdir()), "/")
#     for (ii in 1:length(chrs)) {

#         cat(ii, "/", length(chrs), "...\n", sep = "")

#         xdt <- Rmutmod:::chrom2table(chrs[ii], mafdb, cohort, targetdb, genomePaths[ii], nflank, pkmers, fdirs, fplabs)
#         xdt <- xdt[complete.cases(xdt)]
#         xdt[, ':=' ("start" = NULL, "rangeid" = NULL, "ref" = NULL)]

#         xdt[, "mutcat" := stringi::stri_join(kmer, mut, sep = ">")]
#         xdt[, ':=' ("kmer" = NULL, "mut" = NULL)]
#         Rmutmod:::formatFeatures(xdt, .fplabs)
        
#         if (ii == 1) {
            
#             dsk <- bigReg::data_frame(xdt, 1000000, d, 1)

#         } else {

#             dsk$append(xdt)
            
#         }
#         rm(xdt)

#     }

#     formula_string <- paste0("nmut ~ mutcat * (", paste(names(fplabs), collapse = " + "), ")")
#     model <- bigReg::bglm(
#         formula = as.formula(formula_string),
#         family = bigReg::poisson_(),
#         data = dsk,
#         control = bigReg:::.control()
#     )

#     mutGLMs <- new_MutCGLM(
#         model = model,
#         mafdir = mafdir,
#         cohort = cohort,
#         k = k,
#         targetdir = targetdir,
#         genomedir = genomedir,
#         chrs = chrs,
#         fdirs = fdirs,
#         fplabs = fplabs,
#         formula = as.formula(formula_string),
#     )

#     return(mutGLMs)

# }

