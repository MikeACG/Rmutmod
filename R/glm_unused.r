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

# vector2sdt <- function(x, jcum) {

#     UseMethod("vector2sdt", x)

# }

# vector2sdt.factor <- function(x, jcum) {

#     # initialize data table
#     .sdt <- data.table::data.table("i" = 1:length(x))

#     # add dummy column where a 1 should be placed for each observation
#     .sdt[, "j" := match(x, levels(x))]

#     # remove observations of last dummy column and offset column index
#     .sdt <- .sdt[j != length(levels(x))]
#     .sdt <- .sdt[, "j" := j + jcum]

#     # remaining indices should have a 1 as their values
#     .sdt[, "v" := 1L]

#     return(.sdt)

# }

# vector2sdt.default <- function(x, jcum) {

#     # initialize data table
#     .sdt <- data.table::data.table("i" = 1:length(x))

#     # add column where values will go offseted by previous columns
#     # also add the values themselves
#     .sdt[, ':=' ("j" = jcum + 1L, "v" = x)]

#     # remove observations with value of 0
#     .sdt <- .sdt[v != 0]

#     return(.sdt)

# }

# interaction2sdt <- function(x, jcum, mutcats) {

#     UseMethod("interaction2sdt", x)

# }

# interaction2sdt.factor <- function(x, jcum, mutcats) {

#     # initialize data table
#     .sdt <- data.table::data.table("i" = 1:length(x))

#     # add in what level each observation is for mutcats & remove those in last level
#     .sdt[, "jmutcat" := match(mutcats, levels(mutcats))]

#     # add in what level each observation is for the target feature & remove those in last level
#     .sdt[, "jx" := match(x, levels(x))]

#     # remove observations from last level of any of mutcats or x
#     .sdt <- .sdt[jmutcat != length(levels(mutcats)) & jx != length(levels(x))]

#     # jx - 1: "sets" of ncolsMutcats "behind" the set where the 1 should go
#     # the position inside the set is given by jmutcat
#     ncolsMutcats <- length(levels(mutcats)) - 1
#     .sdt[, "j" := jmutcat + (ncolsMutcats * (jx - 1L))]

#     # remove temporary columns and offset all column indices by the cumulative of columns
#     .sdt[, ':=' ("jmutcat" = NULL, "jx" = NULL)]
#     .sdt <- .sdt[, "j" := j + jcum]

#     # remaining indices should have a 1 as their values
#     .sdt[, "v" := 1L]

#     return(.sdt)

# }

# interaction2sdt.default <- function(x, jcum, mutcats) {

#     # initialize data table
#     .sdt <- data.table::data.table("i" = 1:length(x))

#     # add in what level each observation is for mutcats
#     .sdt[, "j" := match(mutcats, levels(mutcats))]

#     # add value of each observation for the target feature
#     .sdt[, "v" := x]

#     # remove observation if in last level of mutcats or a value of 0 for x
#     .sdt <- .sdt[j != length(levels(mutcats)) & v != 0]

#     # offset all column indices by the cumulative of columns
#     .sdt <- .sdt[, "j" := j + jcum]

#     return(.sdt)

# }

# dt2sdt <- function(xdt, mutcatInteractions = TRUE) {

#     # calculate the number of columns that each feature will create in the model matrix
#     nf <- sapply(xdt, function(x) if (is.factor(x)) length(levels(x)) - 1 else 1)
#     nfc <- head(cumsum(c(0, nf)), -1)

#     # get sparse representation
#     sdt <- mapply(vector2sdt, xdt, nfc, SIMPLIFY = FALSE)
#     sdt <- data.table::rbindlist(sdt)

#     # deal with interactions if neccessary
#     isdt <- data.table::data.table()
#     nfi <- 0
#     if (mutcatInteractions) {

#         mutcats <- xdt$mutcat
#         xdt[, "mutcat" := NULL]
#         nfi <- nf[names(xdt)] * (length(levels(mutcats)) - 1)
#         nfic <- head(cumsum(c(0, nfi)), -1) + sum(nf)

#         isdt <- mapply(interaction2sdt, xdt, nfic, MoreArgs = list(mutcats = mutcats), SIMPLIFY = FALSE)
#         isdt <- data.table::rbindlist(isdt)

#     }

#     return(rbind(sdt, isdt))

# }

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

