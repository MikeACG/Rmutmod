#' @import data.table
#' @importFrom dplyr %>%

chrom2mafDesign <- function(.chr, mafdb, targetdb, genomePath, nflank, fdirs) {

    # get kmers of target sites and mutations
    genome <- setNames(Biostrings::readDNAStringSet(genomePath), .chr)
    xdt <- target2xdt(targetdb, .chr, nflank, genome)
    mutdt <- Rmutmod:::maf2mutdt(mafdb, "all", .chr, nflank, genome, setNames("Cohort", "cohort"))
    rm(genome)

    # ensure the observed mutations are falling in the specified target
    mutdt[xdt, "rangeid" := i.rangeid, on = c("start", "mut")]
    mutdt <- mutdt[!is.na(rangeid)]

    # add the model features to each possible mutaiton in the target and observed mutations
    Rmutmod:::addFeatures(fdirs, .chr, xdt, mutdt)
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
        ccxdt <- split(ccxdt, paste(substr(ccxdt$kmer, icenter, icenter), ccxdt$mut, sep = "_"))
        mapply(
            data.table::fwrite,
            ccxdt,
            paste0(tmpdir, names(ccxdt), ".tmp"),
            MoreArgs = list(sep = "\t", append = TRUE, nThread = 1)
        )
        rm(ccxdt); gc()

    }

    return(tmpdir)

}

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
            rank_check = "skip"
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

varTargetFilter <- function(ccxdt, modvars, minTarget) {

    badIdxs <- c()
    for (v in modvars) {

        countdt <- ccxdt[, list("target" = sum(nchance)), by = v]
        print(countdt)
        badLevels <- countdt[target < minTarget][[v]]
        badIdxs <- append(badIdxs, which(ccxdt[[v]] %in% badLevels))

    }

    if (length(badIdxs) > 0) return(ccxdt[-unique(badIdxs)])
    return(ccxdt)

}

# formatFeatures <- function(ccxdt, vars2format) {

#     novarVars <- c()
#     for (jj in 1:length(vars2format)) {

#         v <- vars2format[jj]
#         if (is.factor(ccxdt[[v]]) | is.character(ccxdt[[v]])) {

#             ccxdt[, (v) := factor(as.character(get(v)))]
#             check <- length(levels(ccxdt[[v]])) > 1

#         } else { # numeric or integer

#             check <- var(ccxdt[[v]]) > 0

#         }

#         if (!check) novarVars <- c(novarVars, v)

#     }

#     return(novarVars)

# }

