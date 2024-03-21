#' @import data.table
#' @importFrom dplyr %>%

new_MutMatrix <- function(
    modeldt = data.table::data.table(
        kmer = character(0L),
        ref = character(0L),
        mut = character(0L),
        n = integer(0L),
        abundance.adj = numeric(0L),
        mutRate = numeric(0L)
    ),
    mafdir = character(1L),
    cohort = character(1L),
    k = integer(1L),
    targetdir = character(1L),
    genomedir = character(1L),
    chrs = character(0L)
) {

    mutMatrix <- structure(
        list(
            modeldt = modeldt,
            mafdir = mafdir,
            chort = cohort,
            k = k,
            targetdir = targetdir,
            genomedir = genomedir,
            chrs = chrs
        ),
        class = c("Rmutmod", "MutMatrix")
    )

    return(mutMatrix)

}

#' @export
kGet <- function(x) {

    UseMethod("kGet")

}

kGet.Rmutmod <- function(rmutmod) {

    return(rmutmod$k)

}

makePkmers <- function(k) {

    kmersComp <- gtools::permutations(4, k, c("A", "C", "G", "T"), repeats.allowed = TRUE)
    kmersComp <- kmersComp[kmersComp[, ((k - 1) / 2) + 1] %in% c("C", "T"), ]
    pkmers <- apply(kmersComp, 1, paste, collapse = "")

    return(pkmers[order(kmersComp[, ceiling(k / 2)])])

}

# checks if central nucleotide of kmer is a purine
isPuri <- function(kmers, nflank) {

    centerIdxs <- rep(nflank + 1, length(kmers))
    centers <- Biostrings::subseq(kmers, centerIdxs, centerIdxs)
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
        position = siteIdxs,
        kmer = as.character(genome[kmerRanges], use.names = FALSE),
        rangeid = rangeids
    )

    kmerdt[
        isPuri(kmer, nflank),
        "kmer" := as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(kmer)), use.names = FALSE)
    ]
    
    return(kmerdt)

}

# count pyrimidine oriented kmers in sites where its possible to detect a mutation
countSites <- function(targetdb, .chr, nflank, genome, pkmers) {

    targetdt <- targetdb %>% 
        dplyr::filter(seqnames == .chr) %>%
        dplyr::select(dplyr::all_of(c("start", "end"))) %>%
        dplyr::collect()
    
    kmerdt <- ranges2kmerdt(targetdt$start, targetdt$end, .chr, nflank, genome)
    
    countdt <- data.table::CJ(kmer = pkmers, unique = TRUE)
    countdt[
        kmerdt[, list("abundance" = .N), by = "kmer"],
        "abundance" := i.abundance,
        on = "kmer"
    ]

    return(countdt)

}

countMutations <- function(mafdb, cohort, .chr, nflank, genome, pkmers) {

    # load mutations and remove flagged records
    mafdt <- mafdb %>% 
        dplyr::filter(Cohort == cohort, Chromosome == .chr) %>%
        dplyr::select(dplyr::all_of(c("Start_Position", "Tumor_Seq_Allele2", "modelExclude"))) %>%
        dplyr::collect()
    mafdt <- mafdt[modelExclude == FALSE]

    # if no mutations return 0 counts for all categories
    countdt <- data.table::CJ(kmer = pkmers, mut = c("A", "C", "G", "T"), unique = TRUE)
    if (nrow(mafdt) == 0L) {

        countdt[, "n" := rep(0L, nrow(countdt))]
        return(countdt)

    }

    mafRanges <- GenomicRanges::GRanges(
        rep(.chr, nrow(mafdt)),
        IRanges::IRanges(mafdt$Start_Position - nflank, mafdt$Start_Position + nflank)
    )
    mutdt <- data.table::data.table(
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

    countdt[
        mutdt[, list("n" = .N), by = c("kmer", "mut")],
        "n" := i.n,
        on = c("kmer", "mut")
    ]
    
    return(countdt)

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
        data.table(gendt)[, list("n" = .N), by = "gender"],
        "n" := i.n,
        on = "gender"
    ]
    countdt[is.na(n), "n" := 0L]

    # we will see the abundance counts as being composed of
    # equal tumor contributions. If all tumors had the same set of chromosomes,
    # each one would contribute 2 chromosomes per chromosome ID
    adjdt <- data.table(
        abundance = abundance,
        chr = .chrs,
        nchr = 2L * sum(countdt$n)
    )
    adjdt[, "chrContrib" := abundance / nchr]

    # in the case of autosomes, all tumors do contribute 2 chromosomes
    # per chromosome id, but for sex chromosomes this must be adjusted
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
trainMutMat <- function(mafdir, cohort, k, targetdir, genomedir, chrs = c(paste0("chr", 1:22), "chrX", "chrY")) {

    pkmers <- makePkmers(k)
    nflank <- (k - 1) / 2
    mafdb <- arrow::open_dataset(mafdir)
    targetdb <- arrow::open_dataset(targetdir)
    genomePaths <- paste0(genomedir, chrs, ".fasta")
    abudt <- list()
    nmutdt <- list()
    for (ii in 1:length(chrs)) {

        cat(ii, "/", length(chrs), "...\n", sep = "")

        # load genome, count abundances & mutations
        genome <- setNames(Biostrings::readDNAStringSet(genomePaths[ii]), chrs[ii])
        abudt[[ii]] <- countSites(targetdb, chrs[ii], nflank, genome, pkmers)
        nmutdt[[ii]] <- countMutations(mafdb, cohort, chrs[ii], nflank, genome, pkmers)

        rm(genome); gc()

    }
    nmutdt <- data.table::rbindlist(nmutdt)

    # drop mutation count entries with same ref and mut (no mutation)
    icenter <- nflank + 1
    nmutdt[, "ref" := substr(kmer, icenter, icenter)]
    nmutdt <- nmutdt[ref != mut]

    # remaining NA entries are non-found mutations, should be 0. Aggregate across chromosomes
    nmutdt[is.na(n), "n" := 0L]
    nmutdt <- nmutdt[, list("n" = sum(n)), by = c("kmer", "ref", "mut")]
    
    # get gender-adjusted abundances aggregated across chromosomes
    abudt <- cbind(data.table::rbindlist(abudt), "chr" = rep(chrs, sapply(abudt, nrow)))
    abudt[, "abundance.adj" := adjustByGender(abundance, chr, mafdb, cohort)]
    abudt <- abudt[, list("abundance.adj" = sum(abundance.adj)), by = "kmer"]

    # add abundances to mutation counts, get mutation rate
    nmutdt[abudt, "abundance.adj" := i.abundance.adj, on = "kmer"]
    nmutdt[, "mutRate" := n / abundance.adj]

    # make mutational matrix object
    mutMatrix <- new_MutMatrix(
        nmutdt,
        mafdir,
        cohort,
        k,
        targetdir,
        genomedir,
        chrs
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

modelGet.MutMatrix <- function(mutmatrix) {

    return(mutmatrix$modeldt)

}

#' @export
mutpredict <- function (x, ...) {

   UseMethod("mutpredict", x)

}

mutpredict.MutMatrix <- function(mutmatrix, newdata, ...) {

    modeldt <- modelGet(mutmatrix)
    newdata[modeldt, "mutRate" := i.mutRate, on = c("kmer", "mut")]

    return()

}

#' @export
pkmersGet <- function(x) {

    UseMethod("pkmersGet", x)

}

pkmersGet.MutMatrix <- function(mutmatrix) {

    modeldt <- modelGet(mutmatrix)
    return(unique(modeldt$kmer))

}
