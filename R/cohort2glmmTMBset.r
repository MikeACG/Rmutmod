#' @export
cohort2xdt <- function(mafdb, cohort, k, targetdb, genomePaths, chrs, fdirs) {

    nflank <- (k - 1) / 2
    xdt <- list()
    for (ii in 1:length(chrs)) {

        cat(ii, "/", length(chrs), "...\n", sep = "")

        xdt[[ii]] <- chrom2table(chrs[ii], mafdb, cohort, targetdb, genomePaths[ii], nflank, fdirs)

    }
    xdt <- rbindlist(xdt)

    xdt[, "mutcat" := stringi::stri_join(kmer, ">", mut)]
    xdt[, ':=' ("kmer" = NULL, "mut" = NULL, "start" = NULL, "rangeid" = NULL)]

    return(xdt)

}

#' @export
cohort2files <- function(mafdb, cohort, k, targetdb, genomePaths, chrs, fdirs) {

    # create directory for temporal files to avoid using a lot of memory
    tmpdir <- paste0("./", basename(tempdir()), "/")
    if (file.exists(tmpdir)) unlink(tmpdir, recursive = TRUE)
    dir.create(tmpdir)

    nflank <- (k - 1) / 2
    icenter <- nflank + 1
    for (ii in 1:length(chrs)) {

        cat(ii, "/", length(chrs), "...\n", sep = "")

        xdt <- chrom2table(chrs[ii], mafdb, cohort, targetdb, genomePaths[ii], nflank, fdirs)
        xdt2disk(xdt, tmpdir)
        rm(xdt)

    }

    return(tmpdir)

}

#' @export
xdt2monoGLMMTMB <- function(xdt, .cond) {

    # format the variables correctly
    for (v in names(xdt)) {

        if (!is.numeric(xdt[[v]])) xdt[, (v) := factor(as.character(get(v)))]

    }

    .REML <- FALSE
    if (any(grepl("[|]", as.character(.cond)))) .REML <- TRUE
    model <- glmmTMB::glmmTMB(
        .cond,
        xdt,
        glmmTMB::nbinom2(),
        dispformula = ~ 1,
        sparseX = c("cond" = TRUE, "zi" = FALSE, "disp" = TRUE),
        control = glmmTMB::glmmTMBControl(
            optCtrl = list(iter.max = 10000, eval.max = 10000),
            #optimizer = optim,
            #optArgs = list(method = "BFGS")
        ),
        REML = FALSE,
        verbose = TRUE
    )

    return(model)

}


