#' @import data.table
#' @importFrom dplyr %>%

#' @export
files4monoGLMMTMB <- function(mafdb, k, targetdb, genomePaths, chrs, fdirs) {

    # create directory for temporal files to avoid using a lot of memory
    tmpdir <- paste0("./", basename(tempdir()), "/")
    if (file.exists(tmpdir)) unlink(tmpdir, recursive = TRUE)
    dir.create(tmpdir)

    nflank <- (k - 1) / 2
    icenter <- nflank + 1
    for (ii in 1:length(chrs)) {

        cat(ii, "/", length(chrs), "...\n", sep = "")

        # get possible mutations in this chromosome and observed mutation counts
        ccxdt <- chrom2mafDesign(chrs[ii], mafdb, targetdb, genomePaths[ii], nflank, fdirs)

        # save to disk by mononucleotide substitution type
        ccxdt <- split(ccxdt, paste(substr(xdt$kmer, icenter, icenter), xdt$mut, sep = "_"))
        mapply(
            data.table::fwrite,
            ccxdt,
            paste0(tmpdir, names(xdt), ".tmp"),
            MoreArgs = list(sep = "\t", append = TRUE, nThread = 1)
        )
        rm(xdt); gc()

    }

    return(tmpdir)

}

#' @export
file2monoGLMMTMB <- function(tmpPath, .cond, ntumordt) {

    # aggregate data by design features
    ccxdt <- data.table::fread(tmpPath, nThread = 1)
    modvars <- setdiff(names(ccxdt), c("nmut", "nchance", "mut"))
    ccxdt <- ccxdt[, list("nmut" = sum(nmut), "nchance" = sum(nchance)), by = modvars] # this does nothing when 1 of the variables doesnt repeat across chromosomes (frequent with random effects)

    # annotate with the number of samples per cohort and adjust abundances by them
    ccxdt[ntumordt, "ntumor" := i.ntumor, on = "cohort"]
    ccxdt[, "nchanceAdj" := nchance * ntumor] # mutation chances
    ccxdt[, "logchance" := log(nchanceAdj)]

    # format the variables correctly and check that they have variance
    ccxdt[, ':=' ("nchance" = NULL, "nchanceAdj" = NULL, "ntumor" = NULL)]
    novarVars <- formatFeatures(ccxdt, modvars)
    .cond <- remove_mterms(.cond, novarVars) # if variables have no variance, remove them from modelformula

    # this fits for T>G mutations: .cond <- as.formula("nmut ~ kmer + tx + rxMCF7 + nucLBBC_5bins + (kmer || cohort) + (kmer || binGenome10mb) + offset(logchance)")
    model <- glmmTMB::glmmTMB(
        .cond,
        ccxdt,
        glmmTMB::nbinom2(),
        dispformula = ~ 1,
        sparseX = c("cond" = TRUE, "zi" = FALSE, "disp" = TRUE),
        control = glmmTMB::glmmTMBControl("optCtrl" = list(iter.max = 10000, eval.max = 10000)),
        REML = TRUE,
        verbose = TRUE
    )
    rm(ccxdt); gc()

    # running VarrCorr can temporarily bloat the memory. Since the fitting model part could require a big
    # machine, we take aadvantage of that and run VarrCorr right now to store the result
    monoMAFglmmTMB <- new_MonoMAFglmmTMB(model, glmmTMB::VarCorr(model))
    return(monoMAFglmmTMB)

}

#' @export
maf2glmmTMBset <- function(mafdir, k, targetdir, genomedir, chrs, fdirs, .cond, outdir) {

    mafdb <- arrow::open_dataset(mafdir)
    targetdb <- arrow::open_dataset(targetdir)
    genomePaths <- paste0(genomedir, chrs, ".fasta")
    tmpdir <- files4monoGLMMTMB(mafdb, k, targetdb, genomePaths, chrs, fdirs)

    ntumordt <- getMAFntumors(mafdb)
    tmpfiles <- list.files(tmpdir, full.names = TRUE)
    bnames <- sub("[.].*", "", basename(tmpfiles))
    if (!file.exists(outdir)) dir.create(outdir)
    outPaths <- setNames(paste0(outdir, bnames, ".rds"), bnames)
    for (ii in 1:length(tmpfiles)) {

        cat(ii, "/", length(tmpfiles), "...\n", sep = "")

        model <- file2monoGLMMTMB(tmpfiles[ii], .cond, ntumordt)
        saveRDS(model, file = outPaths[ii])
        rm(model); gc()

    }

    mafGLMMTMBset <- new_MultiMAFglmmTMB(
        outPaths,
        mafdir,
        k,
        targetdir,
        genomedir,
        chrs,
        fdirs,
        .cond
    )
    unlink(tmpdir, recursive = TRUE)
    return(mafGLMMTMBset)

}

#' @export
getMAFntumors <- function(mafdb) {

    tumorTab <- mafdb %>%
        dplyr::select(dplyr::all_of(c("cohort" = "Cohort", "tumor" = "Tumor_Sample_Barcode"))) %>%
        dplyr::distinct() %>%
        dplyr::collect()
    ntumordt <- data.table::data.table(tumorTab)[, list("ntumor" = .N), by = "cohort"]

    return(ntumordt)

}

formatFeatures <- function(ccxdt, vars2format) {

    novarVars <- c()
    for (jj in 1:length(vars2format)) {

        v <- vars2format[jj]
        if (is.factor(ccxdt[[v]]) | is.character(ccxdt[[v]])) {

            ccxdt[, (v) := factor(as.character(get(v)))]
            check <- length(levels(ccxdt[[v]])) > 1

        } else { # numeric or integer

            check <- var(ccxdt[[v]]) > 0

        }

        if (!check) novarVars <- c(novarVars, v)

    }

    return(novarVars)

}

# https://stackoverflow.com/questions/23381616/r-update-function-how-to-drop-all-the-variables-that-are-related-to-a-pre-speci
# this currently does not support removing from random effects formula part
remove_terms <- function(form, term) {

  fterms <- terms(form)
  fac <- attr(fterms, "factors")
  idx <- which(as.logical(fac[term, ]))
  new_fterms <- drop.terms(fterms, dropx = idx, keep.response = TRUE)

  return(formula(new_fterms))

}

remove_mterms <- function(form, terms) {

    jj <- 1
    nform <- form
    while(jj <= length(terms)) {

        warning(paste0("Removing ", terms[jj], " from model formula as fit will fail if included."))
        nform <- remove_terms(nform, terms[jj])
        jj <- jj + 1

    }

    return(nform)

}
