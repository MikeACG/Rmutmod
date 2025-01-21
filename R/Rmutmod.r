#' @import data.table
#' @import glmnet
#' @importFrom dplyr %>%

makePkmers <- function(k) {

    kmersComp <- gtools::permutations(4, k, c("A", "C", "G", "T"), repeats.allowed = TRUE)
    kmersComp <- kmersComp[kmersComp[, ((k - 1) / 2) + 1] %in% c("C", "T"), ]
    pkmers <- apply(kmersComp, 1, paste, collapse = "")

    return(pkmers[order(kmersComp[, ceiling(k / 2)])])

}

makePcats <- function(pkmers, nflank) {

    catdt <- Rmutmod:::expandMuts(data.table::data.table("kmer" = pkmers), nflank)
    catdt[, "mutcat" := paste(kmer, mut, sep = ">")]

    return(catdt$mutcat)

}

#' @export
mafLoad <- function(mafdb, .cols, .chr = "all", cohort = "all", .vartype = "all", flaggedMuts = "yes", utx = "all") {

    chrQuery <- mafdb %>% dplyr::filter(Chromosome %in% .chr)
    if (.chr[1] == "all") chrQuery <- mafdb

    cohortQuery <- chrQuery %>% dplyr::filter(Cohort %in% cohort)
    if (cohort[1] == "all") cohortQuery <- chrQuery

    vartypeQuery <- cohortQuery %>% dplyr::filter(Variant_Classification %in% .vartype)
    if (.vartype[1] == "all") vartypeQuery <- cohortQuery
    if (.vartype[1] == "exonic") vartypeQuery <- cohortQuery %>% dplyr::filter(!is.na(Variant_Classification))

    excludeQuery <- vartypeQuery
    if (flaggedMuts == "no") excludeQuery <- vartypeQuery %>% dplyr::filter(modelExclude == FALSE) 

    txQuery <- excludeQuery %>% dplyr::filter(Transcript_ID %in% utx)
    if (utx[1] == "all") txQuery <- excludeQuery

    mafdt <- txQuery %>%
        dplyr::select(dplyr::all_of(.cols)) %>%
        dplyr::collect()

    return(mafdt)

}

#' @export
filterTxCount <- function(mafdt, minmut) {

    mafdt[
        mafdt[, list("txmut" = .N), by = "Transcript_ID"],
        "txmut" := i.txmut,
        on = "Transcript_ID"
    ]
    
    return(mafdt[txmut >= minmut])

}

# checks if central nucleotide of kmer is a purine
#' @export
isPuri <- function(kmers, nflank) {

    centerIdxs <- rep(nflank + 1, length(kmers))
    centers <- substr(kmers, centerIdxs, centerIdxs)
    ispuri <- centers == "A" | centers == "G"

    return(ispuri)

}

#' @export
pyriOrient <- function(.dt, .isPuri, icols, ocols) {

    .dt[, (ocols) := .SD, .SDcols = icols]
    .dt[
        .isPuri == TRUE,
        (ocols) := lapply(
            .SD,
            function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
        ),
        .SDcols = icols
    ]

}

#' @export
ranges2kmerdt <- function(.start, .end, .chr, nflank, genome) {

    siteIdxs <- mapply(':', .start, .end, SIMPLIFY = FALSE)
    rangeids <- rep(1:length(siteIdxs), sapply(siteIdxs, length))
    siteIdxs <- unlist(siteIdxs, recursive = FALSE, use.names = FALSE)

    kmerRanges <- GenomicRanges::GRanges(rep(.chr, length(siteIdxs)), IRanges::IRanges(siteIdxs - nflank, siteIdxs + nflank))
    kmerdt <- data.table::data.table(
        start = siteIdxs,
        kmer = as.character(genome[kmerRanges], use.names = FALSE),
        rangeid = rangeids
    )
    
    return(kmerdt)

}

# get pyrimidine oriented kmers in sites where its possible to detect a mutation
target2kmerdt <- function(targetdb, .chr, nflank, genome) {

    targetdt <- targetdb %>% 
        dplyr::filter(seqnames == .chr) %>%
        dplyr::select(dplyr::all_of(c("start", "end"))) %>%
        dplyr::collect()
    
    tkmerdt <- ranges2kmerdt(targetdt$start, targetdt$end, .chr, nflank, genome)
    pyriOrient(tkmerdt, isPuri(tkmerdt$kmer, nflank), "kmer", "kmer")

    return(tkmerdt)

}

maf2mutdt <- function(mafdb, cohort, .chr, nflank, genome, extraCols = c()) {

    # load mutations without flagged records
    cols <- c(setNames(c("Start_Position", "Tumor_Seq_Allele2"), c("start", "mut")), extraCols)
    mafdt <- mafLoad(mafdb, cols, .chr, cohort, flaggedMuts = FALSE)

    # if no mutations return 0 counts for all categories
    mutdt <- data.table::data.table(start = integer(0), kmer = character(0), mut = character(0), Cohort = character(0))
    if (nrow(mafdt) == 0L) return(mutdt)

    mafRanges <- GenomicRanges::GRanges(
        rep(.chr, nrow(mafdt)),
        IRanges::IRanges(mafdt$start - nflank, mafdt$start + nflank)
    )
    mutdt <- mafdt[, .SD, .SDcols = names(cols)]
    mutdt[, "kmer" := as.character(genome[mafRanges], use.names = FALSE)]

    mutdt[
        isPuri(kmer, nflank) == TRUE,
        ':=' (
            "kmer" = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(kmer)), use.names = FALSE),
            "mut" = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(mut)), use.names = FALSE)
        )
    ]
    
    return(mutdt)

}

addFeature <- function(tdt, featuredt, fname) {

    if (nrow(tdt) == 0) {

        tdt[, (fname) := character(0L)]
        return()

    }

    tdt[
        featuredt,
        (fname) := i.feature,
        on = "start"
    ]

}

addFeatures <- function(fdirs, .chr, ...) {

    dtlist <- list(...)

    jj <- 1L
    while (jj <= length(fdirs)) {

        featuredb <- arrow::open_dataset(fdirs[jj])
        featuredt <- featuredb %>% 
            dplyr::filter(seqnames == .chr) %>%
            dplyr::collect()

        lapply(dtlist, addFeature, featuredt, names(fdirs)[jj])

        jj <- jj + 1
        rm(featuredt, featuredb)

    }

}

dtcount <- function(.dt, pvl, cname) {

    countdt <- do.call(
        data.table::CJ,
        c(pvl, unique = TRUE)
    )
    countdt[
        .dt[, list("count" = .N), by = names(pvl)],
        (cname) := i.count,
        on = names(pvl)
    ]

    return(countdt)

}

chrom2matrix <- function(.chr, mafdb, cohort, targetdb, genomePath, nflank, pkmers, fdirs, fplabs) {

    genome <- setNames(Biostrings::readDNAStringSet(genomePath), .chr)
    tkmerdt <- target2kmerdt(targetdb, .chr, nflank, genome)
    mutdt <- maf2mutdt(mafdb, cohort, .chr, nflank, genome)
    rm(genome)

    addFeatures(fdirs, .chr, tkmerdt, mutdt)

    pAbu <- append(list("kmer" = pkmers), fplabs)
    abudt <- dtcount(tkmerdt, pAbu, "abundance")

    pMut <- append(list("mut" = c("A", "C", "G", "T")), pAbu)
    matrixdt <- dtcount(mutdt, pMut, "n")

    matrixdt[abudt, "abundance" := i.abundance, on = names(pAbu)]
    return(matrixdt)

}

adjustByGender <- function(abundance, .chrs, mafdb, cohort) {

    # get gender of unique tumors in cohort
    gendt <- mafdb %>% 
        dplyr::filter(Cohort == cohort) %>%
        dplyr::distinct(Tumor_Sample_Barcode, gender) %>%
        dplyr::select("gender") %>%
        dplyr::collect()
    
    # count gender of tumors
    countdt <- data.table::CJ(gender = c("UKNOWN", "FEMALE", "MALE"))
    countdt[
        data.table::data.table(gendt)[, list("n" = .N), by = "gender"],
        "n" := i.n,
        on = "gender"
    ]
    countdt[is.na(n), "n" := 0L]

    # we will see the abundance counts as being composed of
    # equal tumor contributions. If all tumors had the same set of chromosomes,
    # each one would contribute 2 chromosomes per chromosome ID
    adjdt <- data.table::data.table(
        abundance = abundance,
        chr = .chrs,
        nchr = 2L * sum(countdt$n)
    )
    adjdt[, "chrContrib" := abundance / nchr]

    # in the case of autosomes, all tumors do contribute 2 chromosomes
    # per chromosome id, but for sex chromosomes this must be adjusted
    # for "UKNOWN" gender tumors, it is assumed mutations in the sex chromosomes were excluded
    adjdt[, "nchr.adj" := nchr]
    adjdt[
        chr == "chrX",
        "nchr.adj" := (2L * countdt[gender == "FEMALE"]$n) + countdt[gender == "MALE"]$n
    ]
    adjdt[
        chr == "chrY",
        "nchr.adj" := countdt[gender == "MALE"]$n
    ]
    adjdt[, "abundance.adj" := chrContrib * nchr.adj]

    return(adjdt$abundance.adj)

}

#' @export
trainMutMat <- function(mafdir, cohort, k, targetdir, genomedir, chrs, fdirs, fplabs) {

    pkmers <- makePkmers(k)
    nflank <- (k - 1) / 2
    mafdb <- arrow::open_dataset(mafdir)
    targetdb <- arrow::open_dataset(targetdir)
    genomePaths <- paste0(genomedir, chrs, ".fasta")
    matrixdt <- list()
    for (ii in 1:length(chrs)) {

        cat(ii, "/", length(chrs), "...\n", sep = "")

        matrixdt[[ii]] <- chrom2matrix(chrs[ii], mafdb, cohort, targetdb, genomePaths[ii], nflank, pkmers, fdirs, fplabs)

    }
    nrchr <- sapply(matrixdt, nrow)
    matrixdt <- data.table::rbindlist(matrixdt)
    matrixdt[, "seqnames" := rep(chrs, nrchr)]

    # drop mutation count entries with same ref and mut (no mutation)
    icenter <- nflank + 1
    matrixdt[, "ref" := substr(kmer, icenter, icenter)]
    matrixdt <- matrixdt[ref != mut]

    # remaining NA entries are non-found mutations or abundances, should be 0
    matrixdt[is.na(n), "n" := 0L]
    matrixdt[is.na(abundance), "abundance" := 0L]
    
    # get gender-adjusted abundances aggregated across chromosomes
    matrixdt[, "abundance.adj" := adjustByGender(abundance, seqnames, mafdb, cohort)]
    #matrixdt[, "abundance.adj" := abundance]

    # aggregate across chromosomes, get mutation rate
    matrixdt <- matrixdt[
        ,
        list("n" = sum(n), "abundance.adj" = sum(abundance.adj)),
        by = c("ref", "kmer", "mut", names(fplabs))
    ]
    matrixdt[, "mutRate" := n / abundance.adj]

    # make mutational matrix object
    mutMatrix <- new_MutMatrix(
        matrixdt,
        mafdir,
        cohort,
        k,
        targetdir,
        genomedir,
        chrs,
        fdirs,
        fplabs
    )

    return(mutMatrix)

}

# this does the same as adjustByGender with the "intuitive" method
# was used to check that produces the same result for a few tumors (2 female, 2 male, 1 uknown)
# this "intuitive" method should not be used because it will be numerically unstable for lots of tumors
sumAbundance2 <- function(abundance, tumordt) {

    xi <- which(names(abundance) == "chrX")
    yi <- which(names(abundance) == "chrY")
    res <- list()
    for (i in 1:nrow(tumordt)) {

        a <- abundance
        if (tumordt$gender[i] == "FEMALE") {

            a[[xi]] <- a[[xi]] * 2L
            a[[yi]] <- a[[yi]] * 0L

        }
        if (tumordt$gender[i] == "UKNOWN") {

            a[[xi]] <- a[[xi]] * 0L
            a[[yi]] <- a[[yi]] * 0L

        }
        a[-c(xi, yi)] <- lapply(a[-c(xi, yi)] ,"*" , 2L)

        res[[i]] <- Reduce("+", a)

    }

    return(matrix(Reduce("+", res), ncol = 1))

}

expandMuts <- function(sitedt, nflank) {

    icenter <- nflank + 1
    sitedt[, "ref" := substr(kmer, icenter, icenter)]
    pmuts <- list("C" = c("A", "G", "T"), "T" = c("A", "C", "G"))[sitedt$ref]
    xdt <- sitedt[rep(1:nrow(sitedt), each = 3L)]
    xdt[, "mut" := unlist(pmuts, recursive = FALSE, use.names = FALSE)]

    sitedt[, "ref" := NULL]
    return(xdt)

}

#' @export
mutdesign <- function(rmutmod, rangedt, .chr) {
    
    genomeDir <- genomedirGet(rmutmod)
    k <- kGet(rmutmod)
    fdirs <- fdirsGet(rmutmod)
    fplabs <- fplabsGet(rmutmod)

    nflank <- (k - 1) / 2
    genome <- setNames(Biostrings::readDNAStringSet(paste0(genomeDir, .chr, ".fasta")), .chr)
    sitedt <- ranges2kmerdt(rangedt$start, rangedt$end, .chr, nflank, genome)
    pyriOrient(sitedt, isPuri(sitedt$kmer, nflank), "kmer", "kmer")
    rm(genome)

    addFeatures(fdirs, .chr, sitedt)
    xdt <- expandMuts(sitedt, nflank)
    formatFeatures(xdt, fplabs)

    return(xdt)

}

formatFeatures <- function(xdt, fplabs) {

    # specify levels of variables
    jj <- 1
    while (jj <= length(fplabs)) {

        cname <- names(fplabs)[jj]
        if (length(fplabs[[jj]]) > 1) xdt[, (cname) := factor(as.character(get(cname)), fplabs[[jj]])]
        jj <- jj + 1

    }

    return()

}

chrom2table <- function(.chr, mafdb, cohort, targetdb, genomePath, nflank, pkmers, fdirs, fplabs) {

    # get kmers of target sites and mutations
    genome <- setNames(Biostrings::readDNAStringSet(genomePath), .chr)
    tkmerdt <- target2kmerdt(targetdb, .chr, nflank, genome)
    mutdt <- maf2mutdt(mafdb, cohort, .chr, nflank, genome)
    rm(genome)

    # add the model features to each site in the target and get all possible mutations
    tkmerdt <- tkmerdt[grepl("N", kmer) == FALSE] # remove invalid nucleotides
    addFeatures(fdirs, .chr, tkmerdt)
    xdt <- expandMuts(tkmerdt, nflank)
    rm(tkmerdt)

    # count observed mutations at each site
    xdt[
        mutdt[, list("nmut" = .N), by = c("start", "mut")],
        "nmut" := i.nmut,
        on = c("start", "mut")
    ]
    xdt[is.na(nmut), "nmut" := 0L]

    return(xdt)

}

target2pmuts <- function(.chr, targetdb, genome, nflank, fdirs) {

    # get kmers of target sites
    tkmerdt <- Rmutmod:::target2kmerdt(targetdb, .chr, nflank, genome)

    # add the model features to each site in the target and get all possible mutations
    tkmerdt <- tkmerdt[grepl("N", kmer) == FALSE] # remove invalid nucleotides
    Rmutmod:::addFeatures(fdirs, .chr, tkmerdt)
    xdt <- Rmutmod:::expandMuts(tkmerdt, nflank)
    
    return(xdt)

}

chrom2mafDesign <- function(.chr, mafdb, targetdb, genomePath, nflank, fdirs, fplabs) {

    # get kmers of target sites and mutations
    genome <- setNames(Biostrings::readDNAStringSet(genomePath), .chr)
    tkmerdt <- Rmutmod:::target2kmerdt(targetdb, .chr, nflank, genome)
    mutdt <- maf2mutdt(mafdb, "all", .chr, nflank, genome, setNames("Cohort", "cohort"))
    rm(genome)

    # add the model features to each site in the target and get all possible mutations
    tkmerdt <- tkmerdt[grepl("N", kmer) == FALSE] # remove invalid nucleotides
    Rmutmod:::addFeatures(fdirs, .chr, tkmerdt, mutdt)
    xdt <- Rmutmod:::expandMuts(tkmerdt, nflank)
    rm(tkmerdt)

    # aggregate possible mutations by feature categories and type
    ccols <- c("kmer", "mut", names(fplabs))
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

annotateCohortCount <- function(xdt, mafdb) {

    .cols <- c("cohort" = "Cohort", "tumor" = "Tumor_Sample_Barcode")
    tumorTab <- mafdb %>%
        dplyr::select(dplyr::all_of(.cols)) %>%
        dplyr::distinct() %>%
        dplyr::collect()
    
    xdt[
        data.table::data.table(tumorTab)[, list("ntumor" = .N), by = "cohort"],
        "ntumor" := i.ntumor,
        on = "cohort"
    ]
}

#' @export
trainMutCPR2 <- function(mafdir, k, targetdir, genomedir, chrs, fdirs, fplabs, .formula) {

    pkmers <- Rmutmod:::makePkmers(k)
    nflank <- (k - 1) / 2
    mafdb <- arrow::open_dataset(mafdir)
    targetdb <- arrow::open_dataset(targetdir)
    genomePaths <- paste0(genomedir, chrs, ".fasta")
    
    ccxdt <- list()
    for (ii in 1:length(chrs)) {

        cat(ii, "/", length(chrs), "...\n", sep = "")

        ccxdt[[ii]] <- chrom2mafDesign(chrs[ii], mafdb, targetdb, genomePaths[ii], nflank, fdirs, fplabs)

    }

    # aggregate data by design features
    ccxdt <- data.table::rbindlist(ccxdt)
    ccxdt <- ccxdt[, list("nmut" = sum(nmut), "nchance" = sum(nchance)), by = c(names(fplabs), "kmer", "mut", "cohort")] # this does nothing when 1 of the variables doesnt repeat across chromosomes (frequent with random effects)

    # annotate with the number of samples per cohort and adjust abundances by hem
    annotateCohortCount(ccxdt, mafdb)
    ccxdt[, "nchanceAdj" := nchance / 1000 * ntumor] # mutation chances per thousand chances
    ccxdt[, "logchance" := log(nchanceAdj)]

    # format the variables correctly
    ccxdt[, "mutcat" := stringi::stri_join(kmer, mut, sep = ">")]
    ccxdt[, ':=' ("kmer" = NULL, "mut" = NULL, "nchance" = NULL, "nchanceAdj" = NULL, "ntumor" = NULL)]
    fplabs$binGenome10mbClust <- unique(ccxdt$binGenome10mbClust)
    Rmutmod:::formatFeatures(
        ccxdt,
        append(
            fplabs,
            list(
                "mutcat" = makePcats(pkmers, nflank),
                "cohort" = unique(ccxdt$cohort)
            )
        )
    )

    .disp <- as.formula("~ 1")
    #.cond <- as.formula("nmut ~ mutcat + nucLBBC_5bins + rxMCF7 + meNEU_5bins + tx + offset(logchance)")
    #.cond <- as.formula("nmut ~ (mutcat * cohort) + tx + rxMCF7 + meNEU_5bins + nucLBBC_5bins + offset(logchance)")
    .cond <- as.formula("nmut ~ (mutcat * cohort) + (mutcat | binGenome10mbClust) + tx + rxMCF7 + meNEU_5bins + nucLBBC_5bins + offset(logchance)")
    model <- glmmTMB::glmmTMB(
        .cond,
        ccxdt,
        glmmTMB::nbinom2(),
        dispformula = .disp,
        sparseX = c("cond" = TRUE, "zi" = FALSE, "disp" = TRUE),
        control = glmmTMB::glmmTMBControl("optCtrl" = list(iter.max = 10000, eval.max = 10000)),
        verbose = TRUE
    )

    mutGLMMTMB <- new_MutGLMMTMB(
        model = model,
        mafdir = mafdir,
        cohort = cohort,
        k = k,
        targetdir = targetdir,
        genomedir = genomedir,
        chrs = chrs,
        fdirs = fdirs,
        fplabs = fplabs,
        formula = .formula
    )

    return(mutGLMMTMB)

}


chrom2ctable <- function(.chr, mafdb, cohort, targetdb, genomePath, nflank, pkmers, fdirs, fplabs) {

    # get kmers of target sites and mutations
    genome <- setNames(Biostrings::readDNAStringSet(genomePath), .chr)
    tkmerdt <- Rmutmod:::target2kmerdt(targetdb, .chr, nflank, genome)
    mutdt <- maf2mutdt(mafdb, cohort, .chr, nflank, genome)
    rm(genome)

    # add the model features to each site in the target and get all possible mutations
    tkmerdt <- tkmerdt[grepl("N", kmer) == FALSE] # remove invalid nucleotides
    Rmutmod:::addFeatures(fdirs, .chr, tkmerdt, mutdt)
    xdt <- Rmutmod:::expandMuts(tkmerdt, nflank)
    rm(tkmerdt)

    # aggregate mutations by combinations of features and combine with cancer type
    ccols <- c("kmer", "mut", names(fplabs))
    cxdt <- xdt[, list("nchance" = .N), by = ccols]
    rm(xdt)

    # count observed mutations at each site
    cxdt[
        mutdt[, list("nmut" = .N), by = ccols],
        "nmut" := i.nmut,
        on = ccols
    ]
    cxdt[is.na(nmut), "nmut" := 0L]

    return(cxdt)

}

xdt2disk <- function(xdt, tmpdir) {

    xlist <- split(xdt, by = c("kmer", "mut"))
    fnames <- paste0(tmpdir, gsub("[.]", "_", names(xlist)), ".tmp")
    .append <- file.exists(fnames)
    for (jj in 1:length(xlist)) {

        .dt <- xlist[[jj]]
        .dt[, ':=' ("start" = NULL, "rangeid" = NULL, "ref" = NULL, "mut" = NULL, "kmer" = NULL)]
        data.table::fwrite(.dt, fnames[jj], append = .append[jj], nThread = 1)

    }

    return()

}

tryglm <- function(xdt, .formula, .trace = FALSE) {

    # first try to fit a negative binomial model
    m <- tryCatch(
        {

            withCallingHandlers(
                MASS::glm.nb(.formula, xdt, control = glm.control(maxit = 100, trace = .trace)),
                warning = function(w) {
                    lastWarn <<- w
                    invokeRestart("muffleWarning")
                }
            )

        },
        error = function(cnd) list()
    )

    # if fit failed, try a poisson model
    if (length(m) == 0) {

        m <- withCallingHandlers(
            glm(.formula, xdt, family = "poisson", control = glm.control(maxit = 100, trace = .trace)),
            warning = function(w) {
                lastWarn <<- w
                invokeRestart("muffleWarning")
            }
        )

    }

    return(m)

}

strip_glm <- function(m1) {

    m1$data <- NULL
    m1$y <- NULL
    m1$linear.predictors <- NULL
    m1$weights <- NULL
    m1$fitted.values <- NULL
    m1$model <- NULL
    m1$prior.weights <- NULL
    m1$residuals <- NULL
    m1$effects <- NULL
    m1$qr$qr <- NULL
    return(m1)

}

#' @export
trainMutGLMs <- function(mafdir, cohort, k, targetdir, genomedir, chrs, fdirs, fplabs, .formula) {

    pkmers <- Rmutmod:::makePkmers(k)
    nflank <- (k - 1) / 2
    mafdb <- arrow::open_dataset(mafdir)
    targetdb <- arrow::open_dataset(targetdir)
    genomePaths <- paste0(genomedir, chrs, ".fasta")
    
    tmpdir <- paste0("./", basename(tempdir()), "/")
    dir.create(tmpdir)
    for (ii in 1:length(chrs)) {

        cat(ii, "/", length(chrs), "...\n", sep = "")

        xdt <- Rmutmod:::chrom2table(chrs[ii], mafdb, cohort, targetdb, genomePaths[ii], nflank, pkmers, fdirs, fplabs)
        Rmutmod:::xdt2disk(xdt, tmpdir)
        rm(xdt)

    }
    
    # fit models
    fnames <- list.files(tmpdir, full.names = TRUE)
    mnames <- gsub("[.]tmp", "", basename(fnames))
    models <- list()
    warns <- setNames(character(length(fnames)), mnames)
    for (ii in 1:length(fnames)) {

        cat(ii, "/", length(fnames), "...\n", sep = "")

        xdt <- data.table::fread(fnames[ii])
        formatFeatures(xdt, fplabs)
        lastWarn <- list(message = "none")
        models[[mnames[ii]]] <- strip_glm(tryglm(xdt, .formula))
        warns[ii] <- lastWarn$message
        rm(xdt)

    }
    unlink(tmpdir, recursive = TRUE)

    mutGLMs <- new_MutGLMs(
        models = models,
        mafdir = mafdir,
        cohort = cohort,
        k = k,
        targetdir = targetdir,
        genomedir = genomedir,
        chrs = chrs,
        fdirs = fdirs,
        fplabs = fplabs,
        formula = .formula,
        warns = warns
    )

    return(mutGLMs)

}

#' @export
mutpredict.MutGLMs <- function(rmutmod, newdata, ...) {

    models <- modelGet(rmutmod)

    mdt <- data.table::data.table(
        "kmer" = gsub("_.*", "", names(models)),
        "mut" = gsub(".*_", "", names(models))
    )
    mdt[, "midx" := 1:nrow(mdt)]
    newdata[mdt, "midx" := i.midx, on = c("kmer", "mut")]
    
    newdata[
        ,
        "mutRate" := predict(
            models[[.BY$midx]],
            .SD,
            "response"
        ),
        by = "midx" 
    ]

    newdata[, "midx" := NULL]
    return()

}

#' @export
rgammaScaler <- function(x, .n) {

   UseMethod("rgammaScaler", x)

}

#' @export
rgammaScaler.negbin <- function(m, .n) {

    .theta <- m$theta
   return(rgamma(.n, .theta) / .theta)

}

#' @export
rgammaScaler.default <- function(m, .n) {

   return(rep(1L, .n))

}

#' @export
trainMutCPR <- function(mafdir, cohort, k, targetdir, genomedir, chrs, fdirs, fplabs, .formula) {

    pkmers <- Rmutmod:::makePkmers(k)
    nflank <- (k - 1) / 2
    mafdb <- arrow::open_dataset(mafdir)
    targetdb <- arrow::open_dataset(targetdir)
    genomePaths <- paste0(genomedir, chrs, ".fasta")
    
    cxdt <- list()
    for (ii in 1:length(chrs)) {

        cat(ii, "/", length(chrs), "...\n", sep = "")

        cxdt[[ii]] <- chrom2ctable(chrs[ii], mafdb, cohort, targetdb, genomePaths[ii], nflank, pkmers, fdirs, fplabs)

    }
    cxdt <- data.table::rbindlist(cxdt)[, list("nmut" = sum(nmut), "nchance" = sum(nchance)), by = c(names(fplabs), "kmer", "mut")]
    cxdt[, "logchance" := log(nchance)]

    # format the variables correctly
    cxdt[, "mutcat" := stringi::stri_join(kmer, mut, sep = ">")]
    cxdt[, ':=' ("kmer" = NULL, "mut" = NULL, "nchance" = NULL)]
    Rmutmod:::formatFeatures(cxdt, append(fplabs, list("mutcat" = makePcats(pkmers, nflank))))
    #contrasts(cxdt$mutcat) <- memisc::contr.sum(catdt$mutcat)

    .disp <- as.formula("~ mutcat + offset(logchance)")
    #.cond <- as.formula("nmut ~ mutcat + nucLBBC_5bins + rxMCF7 + meNEU_5bins + tx + offset(logchance)")
    .cond <- as.formula("nmut ~ mutcat + offset(logchance)")
    model <- glmmTMB::glmmTMB(
        .cond,
        cxdt,
        glmmTMB::nbinom2(),
        dispformula = .disp,
        sparseX = c("cond" = TRUE, "zi" = FALSE, "disp" = TRUE),
        control = glmmTMB::glmmTMBControl("optCtrl" = list(iter.max = 10000, eval.max = 10000)),
        verbose = TRUE
    )

    mutGLMMTMB <- new_MutGLMMTMB(
        model = model,
        mafdir = mafdir,
        cohort = cohort,
        k = k,
        targetdir = targetdir,
        genomedir = genomedir,
        chrs = chrs,
        fdirs = fdirs,
        fplabs = fplabs,
        formula = .formula
    )

    return(mutGLMMTMB)

}

#' @export
mutpredict.MutGLMMTMB <- function(rmutmod, newdata, ...) {

    model <- modelGet(rmutmod)
    blist <- lapply(glmmTMB::fixef(model), function(b) ifelse(is.na(b), 0.0, b))
    .formula <- as.formula(paste("~", as.character(formulaGet(rmutmod))[3]))

    catMissing <- !("mutcat" %in% names(newdata))
    if (catMissing) newdata[, "mutcat" := stringi::stri_join(kmer, mut, sep = ">")]
    offsetMissing <- !("logchance" %in% names(newdata))
    if (offsetMissing) newdata[, "logchance" := 0L]

    k <- kGet(rmutmod)
    fplabs <- fplabsGet(rmutmod)
    formatFeatures(newdata, append(fplabs, list("mutcat" = makePcats(Rmutmod:::makePkmers(k), floor(k / 2)))))
    X <- Matrix::sparse.model.matrix(.formula, newdata)
    newdata[, "mutRate" := model$modelInfo$family$linkinv((X %*% blist$cond)[, 1] + logchance)]

    if (catMissing) newdata[, "mutcat" := NULL]
    if (offsetMissing) newdata[, "logchance" := NULL]
    return()

}

#' @export
trainMutLRGLM <- function(mafdir, k, targetdir, genomedir, chrs, fdirs, fplabs, .formula, regionVar) {

    pkmers <- Rmutmod:::makePkmers(k)
    nflank <- (k - 1) / 2
    mafdb <- arrow::open_dataset(mafdir)
    targetdb <- arrow::open_dataset(targetdir)
    genomePaths <- paste0(genomedir, chrs, ".fasta")
    
    pars <- list()
    for (ii in 1:length(chrs)) {

        cat(ii, "/", length(chrs), "...\n", sep = "")

        ccxdt <- chrom2mafDesign(chrs[ii], mafdb, targetdb, genomePaths[ii], nflank, fdirs, fplabs)

        # annotate with the number of samples per cohort and adjust abundances by them
        annotateCohortCount(ccxdt, mafdb) # getting the source dataframe of sample counts can be done beforehand, not everytime
        ccxdt[, "nchanceAdj" := nchance / 1000 * ntumor] # mutation chances per thousand chances
        ccxdt[, "logchance" := log(nchanceAdj)]

        # format the variables correctly
        ccxdt[, "mutcat" := stringi::stri_join(kmer, mut, sep = ">")]
        ccxdt[, ':=' ("kmer" = NULL, "mut" = NULL, "nchance" = NULL, "nchanceAdj" = NULL, "ntumor" = NULL)]
        formatFeatures(
            ccxdt,
            append(
                fplabs,
                list(
                    "mutcat" = makePcats(pkmers, nflank),
                    "cohort" = unique(ccxdt$cohort)
                )
            )
        )

        # fit parameters
        ccxdt <- split(ccxdt, ccxdt[[regionVar]])
        pars[[ii]] <- lapply(ccxdt, GLMMTMBpars, .cond, .fam)

    }

    
    #contrasts(cxdt$mutcat) <- memisc::contr.sum(catdt$mutcat)

    .disp <- as.formula("~ mutcat + offset(logchance)")
    #.cond <- as.formula("nmut ~ mutcat + nucLBBC_5bins + rxMCF7 + meNEU_5bins + tx + offset(logchance)")

    mutGLMMTMB <- new_MutGLMMTMB(
        model = model,
        mafdir = mafdir,
        cohort = cohort,
        k = k,
        targetdir = targetdir,
        genomedir = genomedir,
        chrs = chrs,
        fdirs = fdirs,
        fplabs = fplabs,
        formula = .formula
    )

    return(mutGLMMTMB)

}

GLMMTMBpars <- function(.xdt, .cond, .fam) {

    model <- glmmTMB::glmmTMB(
        .cond,
        .xdt,
        .fam,
        sparseX = c("cond" = TRUE, "zi" = FALSE, "disp" = FALSE),
        control = glmmTMB::glmmTMBControl("optCtrl" = list(iter.max = 10000, eval.max = 10000)),
        verbose = TRUE
    )

}
