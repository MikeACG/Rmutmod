makePkmers <- function(k) {

    kmersComp <- gtools::permutations(4, k, c("A", "C", "G", "T"), repeats.allowed = TRUE)
    kmersComp <- kmersComp[kmersComp[, ((k - 1) / 2) + 1] %in% c("C", "T"), ]
    pkmers <- apply(kmersComp, 1, paste, collapse = "")

    return(pkmers[order(kmersComp[, ceiling(k / 2)])])

}

makePcats <- function(pkmers, nflank) {

    catdt <- expandMuts(data.table::data.table("kmer" = pkmers), nflank)
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

    centerIdxs <- rep(nflank + 1L, length(kmers))
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

    return()

}

maf2mutdt <- function(mafdb, cohort, .chr, nflank, genome, extraCols = c()) {

    # load mutations without flagged records
    cols <- c(setNames(c("Start_Position", "Tumor_Seq_Allele2"), c("start", "mut")), extraCols)
    mafdt <- mafLoad(mafdb, cols, .chr, cohort, flaggedMuts = "no")

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
mutdesign <- function(x, rangedt, .chr) {

    UseMethod("mutdesign", x)

}

#' @export
mutdesign.default <- function(rmutmod, rangedt, .chr) {

    xdt <- mutdesignBase(rmutmod, rangedt, .chr)

    return(xdt)

}

#' @export
mutdesign.MonoMAFglmmTMB <- function(rmutmod, rangedt, .chr) {

    xdt <- mutdesignBase(rmutmod, rangedt, .chr)
    xdt[, "mutcat" := stringi::stri_join(kmer, mut, sep = ">")]

    return(xdt)

}

mutdesignBase <- function(rmutmod, rangedt, .chr) {

    genomeDir <- genomedirGet(rmutmod)
    k <- kGet(rmutmod)
    fdirs <- fdirsGet(rmutmod)
    nflank <- (k - 1) / 2

    genome <- setNames(Biostrings::readDNAStringSet(paste0(genomeDir, .chr, ".fasta")), .chr)
    sitedt <- ranges2kmerdt(rangedt$start, rangedt$end, .chr, nflank, genome)
    pyriOrient(sitedt, isPuri(sitedt$kmer, nflank), "kmer", "kmer")
    rm(genome)

    addFeatures(fdirs, .chr, sitedt)
    xdt <- expandMuts(sitedt, nflank)

    return(xdt)

}

chrom2table <- function(.chr, mafdb, cohort, targetdb, genomePath, nflank, fdirs) {

    # get kmers of target sites and mutations
    genome <- setNames(Biostrings::readDNAStringSet(genomePath), .chr)
    tkmerdt <- target2kmerdt(targetdb, .chr, nflank, genome)
    mutdt <- maf2mutdt(mafdb, cohort, .chr, nflank, genome)
    rm(genome)

    # add the model features to each site in the target and get all possible mutations
    tkmerdt <- tkmerdt[grepl("N", kmer) == FALSE] # remove invalid nucleotides
    addFeatures(fdirs, .chr, tkmerdt)
    tkmerdt <- tkmerdt[complete.cases(tkmerdt)]
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

chrom2ctable <- function(.chr, mafdb, cohort, targetdb, genomePath, nflank, pkmers, fdirs, fplabs) {

    # get kmers of target sites and mutations
    genome <- setNames(Biostrings::readDNAStringSet(genomePath), .chr)
    tkmerdt <- target2kmerdt(targetdb, .chr, nflank, genome)
    mutdt <- maf2mutdt(mafdb, cohort, .chr, nflank, genome)
    rm(genome)

    # add the model features to each site in the target and get all possible mutations
    tkmerdt <- tkmerdt[grepl("N", kmer) == FALSE] # remove invalid nucleotides
    addFeatures(fdirs, .chr, tkmerdt, mutdt)
    xdt <- expandMuts(tkmerdt, nflank)
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
