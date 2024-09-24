#' @import data.table
#' @importFrom dplyr %>%

new_MutMatrix <- function(
    modeldt = data.table::data.table(),
    mafdir = character(1L),
    cohort = character(1L),
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

#' @export
kGet <- function(x) {

    UseMethod("kGet")

}

#' @export
kGet.Rmutmod <- function(rmutmod) {

    return(rmutmod$k)

}

makePkmers <- function(k) {

    kmersComp <- gtools::permutations(4, k, c("A", "C", "G", "T"), repeats.allowed = TRUE)
    kmersComp <- kmersComp[kmersComp[, ((k - 1) / 2) + 1] %in% c("C", "T"), ]
    pkmers <- apply(kmersComp, 1, paste, collapse = "")

    return(pkmers[order(kmersComp[, ceiling(k / 2)])])

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

# checks if central nucleotide of kmer is a purine
#' @export
isPuri <- function(kmers, nflank) {

    centerIdxs <- rep(nflank + 1, length(kmers))
    centers <- substr(kmers, centerIdxs, centerIdxs)
    ispuri <- centers == "A" | centers == "G"

    return(ispuri)

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

    kmerdt[
        isPuri(kmer, nflank),
        "kmer" := as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(kmer)), use.names = FALSE)
    ]
    
    return(kmerdt)

}

# get pyrimidine oriented kmers in sites where its possible to detect a mutation
target2kmerdt <- function(targetdb, .chr, nflank, genome) {

    targetdt <- targetdb %>% 
        dplyr::filter(seqnames == .chr) %>%
        dplyr::select(dplyr::all_of(c("start", "end"))) %>%
        dplyr::collect()
    
    tkmerdt <- ranges2kmerdt(targetdt$start, targetdt$end, .chr, nflank, genome)

    return(tkmerdt)

}

maf2mutdt <- function(mafdb, cohort, .chr, nflank, genome) {

    # load mutations and remove flagged records
    mafdt <- mafdb %>% 
        dplyr::filter(Cohort == cohort, Chromosome == .chr) %>%
        dplyr::select(dplyr::all_of(c("Start_Position", "Tumor_Seq_Allele2", "modelExclude"))) %>%
        dplyr::collect()
    mafdt <- mafdt[modelExclude == FALSE]

    # if no mutations return 0 counts for all categories
    mutdt <- data.table::data.table(start = integer(0), kmer = character(0), mut = character(0))
    if (nrow(mafdt) == 0L) return(mutdt)

    mafRanges <- GenomicRanges::GRanges(
        rep(.chr, nrow(mafdt)),
        IRanges::IRanges(mafdt$Start_Position - nflank, mafdt$Start_Position + nflank)
    )
    mutdt <- data.table::data.table(
        start = mafdt$Start_Position,
        kmer = as.character(genome[mafRanges], use.names = FALSE),
        mut = mafdt$Tumor_Seq_Allele2
    )

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

#' @export
modelGet <- function(x) {

    UseMethod("modelGet", x)

}

#' @export
modelGet.MutMatrix <- function(mutmatrix) {

    return(mutmatrix$modeldt)

}

#' @export
mutpredict <- function (x, ...) {

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
mutdesign <- function (x, ...) {

   UseMethod("mutdesign", x)

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
mutdesign.MutMatrix <- function(mutmatrix, rangedt, .chr, ...) {
    
    genomeDir <- genomedirGet(mutmatrix)
    k <- kGet(mutmatrix)
    fdirs <- fdirsGet(mutmatrix)

    nflank <- (k - 1) / 2
    genome <- setNames(Biostrings::readDNAStringSet(paste0(genomeDir, .chr, ".fasta")), .chr)
    sitedt <- ranges2kmerdt(rangedt$start, rangedt$end, .chr, nflank, genome)
    rm(genome)

    addFeatures(fdirs, .chr, sitedt)

    xdt <- expandMuts(sitedt, nflank)

    return(xdt)

}

chrom2table <- function(.chr, mafdb, cohort, targetdb, genomePath, nflank, pkmers, fdirs, fplabs) {

    # get kmers of target sites and mutations
    genome <- setNames(Biostrings::readDNAStringSet(genomePath), .chr)
    tkmerdt <- Rmutmod:::target2kmerdt(targetdb, .chr, nflank, genome)
    mutdt <- Rmutmod:::maf2mutdt(mafdb, cohort, .chr, nflank, genome)
    rm(genome)

    # add the model features to each site in the target and get all possible mutations
    tkmerdt <- tkmerdt[grepl("N", kmer) == FALSE] # remove invalid nucleotides
    Rmutmod:::addFeatures(fdirs, .chr, tkmerdt)
    xdt <- expandMuts(tkmerdt, nflank)
    rm(tkmerdt)

    # count observed mutations at each site
    xdt[
        mutdt[, list("nmut" = .N), by = c("start", "mut")],
        "nmut" := i.nmut,
        on = c("start", "mut")
    ]
    xdt[is.na(nmut), "nmut" := 0L]
    rm(mutdt)

    # format mutation category and leave only relevant columns
    xdt[, ':=' ("start" = NULL, "rangeid" = NULL, "ref" = NULL)]
    xdt[, "mutcat" := stringi::stri_join(kmer, mut, sep = ">")]
    xdt[, ':=' ("kmer" = NULL, "mut" = NULL)]

    # specify levels of variables
    catdt <- expandMuts(data.table::data.table(kmer = pkmers), nflank)
    catdt[, "mutcat" := stringi::stri_join(kmer, mut, sep = ">")]
    fpl <- append(list("mutcat" = catdt$mutcat), fplabs)
    for (jj in 1:length(fpl)) {

        cname <- names(fpl)[jj]
        if (length(fpl[[jj]]) > 0) xdt[, (cname) := factor(get(cname), fpl[[jj]])]

    }

    return(xdt)

}

dt2sdisk <- function(xdt, .formula, xfile, yfile) {

    X <- MatrixModels::model.Matrix(.formula, xdt, sparse = TRUE)
    sdt <- data.table::as.data.table(selectMethod("summary", "dgCMatrix")(X))
    .ncol <- ncol(X)
    rm(X)
    data.table::fwrite(sdt, xfile, append = TRUE, nThread = 1, col.names = FALSE, sep = " ")
    nval <- nrow(sdt)
    rm(sdt)

    # save the outcome variable to its own file and remove it from the table
    ydt <- xdt[, .SD, .SDcols = "nmut"]
    data.table::fwrite(ydt, yfile, append = TRUE, nThread = 1, col.names = FALSE)

    # return matrix market header info
    mminfo <- c("nrow" = nrow(xdt), "ncol" = .ncol, "nval" = nval)
    return(mminfo)

}

#' @export
trainMutGLM <- function(mafdir, cohort, k, targetdir, genomedir, chrs, fdirs, fplabs, .formula) {

    pkmers <- makePkmers(k)
    nflank <- (k - 1) / 2
    mafdb <- arrow::open_dataset(mafdir)
    targetdb <- arrow::open_dataset(targetdir)
    genomePaths <- paste0(genomedir, chrs, ".fasta")
    
    xfile <- tempfile()
    yfile <- tempfile()
    mminfos <- list()
    for (ii in 1:length(chrs)) {

        cat(ii, "/", length(chrs), "...\n", sep = "")

        xdt <- chrom2table(chrs[ii], mafdb, cohort, targetdb, genomePaths[ii], nflank, pkmers, fdirs, fplabs)
        mminfos[[ii]] <- dt2sdisk(xdt, .formula, xfile, yfile)

        rm(xdt)

    }
    
    # handle the header of matrix market file
    .nrow <- sum(sapply(mminfos, '[', 1))
    .ncol <- mminfos[[1]][2]
    .nval <- sum(sapply(mminfos, '[', 3))
    h <- paste0(
        "%%MatrixMarket matrix coordinate real general\n",
        paste(.nrow, .ncol, .nval, sep = " ")
    )
    sfile <- tempfile()
    .cmd <- paste0(
        "echo '",
        h,
        "' > ",
        sfile,
        " && cat ",
        xfile,
        " >> ",
        sfile
    )
    system(.cmd)

    X <- Matrix::readMM(sfile)
    y <- data.table::fread(yfile)[[1]]
    glmMut <- glmnet::glmnet(
        X,
        y,
        "poisson",
        lambda = 0,
        standardize = FALSE
    )

    return()

}

