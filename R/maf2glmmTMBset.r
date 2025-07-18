#' @import data.table
#' @importFrom dplyr %>%

#' @export
chrom2mafDesign <- function(.chr, mafdb, targetdb, genomePath, nflank, fdirs) {

    # get kmers of target sites and mutations
    genome <- setNames(Biostrings::readDNAStringSet(genomePath), .chr)
    xdt <- target2xdt(targetdb, .chr, nflank, genome)
    mutdt <- maf2mutdt(mafdb, "all", .chr, nflank, genome, setNames("Cohort", "cohort"))
    rm(genome)

    # ensure the observed mutations are falling in the specified target
    mutdt[xdt, "rangeid" := i.rangeid, on = c("start", "mut")]
    mutdt <- mutdt[!is.na(rangeid)]

    # add the model features to each possible mutaiton in the target and observed mutations
    addFeatures(fdirs, .chr, xdt, mutdt)
    xdt <- xdt[complete.cases(xdt)]
    mutdt <- mutdt[complete.cases(mutdt)]

    # aggregate possible mutations by feature categories and type
    ccols <- c("kmer", "mut", names(fdirs))
    cxdt <- xdt[, list("nchance" = .N), by = ccols]
    rm(xdt)

    # aggregate observed mutations by feature categories and type
    mccols <- c(ccols, "cohort")
    cmutdt <- mutdt[, list("nmut" = .N), by = mccols]
    rm(mutdt)

    # consider possible mutations in each cohort
    ucohorts <- unique(cmutdt$cohort)
    ccxdt <- cxdt[rep(1:nrow(cxdt), each = length(ucohorts))]
    ccxdt[, "cohort" := rep(ucohorts, nrow(cxdt))]
    rm(cxdt)

    # note number of observed mutations in the possible mutations table
    ccxdt[
        cmutdt,
        "nmut" := i.nmut,
        on = mccols
    ]
    ccxdt[is.na(nmut), "nmut" := 0L]

    return(ccxdt)

}

#' @export
maf2mafDesign <- function(mafdb, k, targetdb, genomePaths, chrs, fdirs) {
    
    nflank <- (k - 1) / 2
    icenter <- nflank + 1
    ccxdt <- list()
    for (ii in 1:length(chrs)) {

        cat(ii, "/", length(chrs), "...\n", sep = "")

        # get possible mutations in this chromosome and observed mutation counts
        ccxdt[[ii]] <- Rmutmod:::chrom2mafDesign(chrs[ii], mafdb, targetdb, genomePaths[ii], nflank, fdirs)

    }

    return(ccxdt)

}

# .df <- expand.grid(ref = c("C", "T"), mut = c("A", "C", "G", "T"), stringsAsFactors = F)
# .df <- .df[.df$ref != .df$mut, ]
# ccxdt <- list()
# for (ii in 1:nrow(.df)) {

#     x <- data.table::fread(paste0("RtmpisxKsy/", .df$ref[ii], "_", .df$mut[ii], ".tmp"), nThread = 1)
#     x[, "mutcat" := stringi::stri_join(kmer, ">", mut)]
#     ccxdt[[ii]] <- x[, list("nmut" = sum(nmut), "nchance" = sum(nchance)), by = setdiff(names(x), c("nmut", "nchance", "mut", "kmer"))]
#     rm(x)

# }
# ccxdt <- data.table::rbindlist(ccxdt)
# modvars <- setdiff(names(ccxdt), c("nmut", "nchance"))
# .cond <- as.formula("nmut ~ mutcat + txSimply + rxMCF7 + rightFlankNuc2 + rightFlankNuc3 + leftFlankNuc2 + leftFlankNuc3 + nucLBBC_5bins + meNEU_5bins + (mutcat + txSimply + rxMCF7 + rightFlankNuc2 + rightFlankNuc3 + leftFlankNuc2 + leftFlankNuc3 + nucLBBC_5bins || cohort) + offset(logchance)")

#' @export
file2monoGLMMTMB <- function(tmpPath, .cond, ntumordt, minTarget = NA_integer_) {

    # aggregate data by design features
    ccxdt <- data.table::fread(tmpPath, nThread = 1)
    modvars <- setdiff(names(ccxdt), c("nmut", "nchance", "mut"))
    ccxdt <- ccxdt[, list("nmut" = sum(nmut), "nchance" = sum(nchance)), by = modvars] # this does nothing when 1 of the variables doesnt repeat across chromosomes (frequent with random effects)

    # annotate with the number of samples per cohort and adjust abundances by them
    ccxdt[ntumordt, "ntumor" := i.ntumor, on = "cohort"]
    ccxdt[, "nchanceAdj" := nchance * ntumor] # mutation chances
    ccxdt[, "logchance" := log(nchanceAdj)]

    # format the variables correctly and check that they have variance
    if (!is.na(minTarget)) ccxdt <- varTargetFilter(ccxdt, modvars, minTarget)
    ccxdt[, ':=' ("nchance" = NULL, "nchanceAdj" = NULL, "ntumor" = NULL)]
    ccxdt[, (modvars) := lapply(.SD, function(x) factor(as.character(x))), .SDcols = modvars]

    .REML <- FALSE
    if (any(grepl("[|]", as.character(.cond)))) .REML <- TRUE
    model <- glmmTMB::glmmTMB(
        .cond,
        ccxdt,
        glmmTMB::nbinom2(),
        dispformula = ~ 1,
        sparseX = c("cond" = TRUE, "zi" = FALSE, "disp" = TRUE),
        control = glmmTMB::glmmTMBControl(
            "optCtrl" = list(iter.max = 10000, eval.max = 10000),
            #rank_check = "skip"
        ),
        REML = .REML,
        verbose = TRUE
    )

    return(model)

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

    multiMAFglmmTMB <- new_MultiMAFglmmTMBparts(
        outPaths,
        mafdir,
        k,
        targetdir,
        genomedir,
        chrs,
        fdirs
    )
    unlink(tmpdir, recursive = TRUE)
    return(multiMAFglmmTMB)

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


