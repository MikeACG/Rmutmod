# this version works in conjuctions with th Rmutmod package

# the idea is to simulate mutations in a "single" genome but in reality they will be the sum of mutations
# in a complete cancer cohort. For this we need the target regions file of sequencing, the genome,
# and the blueprint of the mutational model to use for simulation. The model object should already have all of this

#' @export
mafsim <- function(rmutmod) {

    # get necessary info to simulate from this model
    genomeDir <- genomedirGet(rmutmod)
    k <- kGet(rmutmod)
    fdirs <- fdirsGet(rmutmod)
    fplabs <- fplabsGet(rmutmod)
    targetdir <- targetdirGet(rmutmod)
    chrs <- chrsGet(rmutmod)

    # open neccessary databases and loop the chromosomes
    targetdb <- arrow::open_dataset(targetdir)
    mafdt <- list()
    for (ii in 1:length(chrs)) {

        cat(ii, "/", length(chrs), "...\n", sep = "")

        # get the model frame for each possible mutation
        targetdt <- targetdb %>% 
            dplyr::filter(seqnames == chrs[ii]) %>%
            dplyr::select(dplyr::all_of(c("start", "end"))) %>%
            dplyr::collect()
        xdt <- mutFrame(targetdt, k, genomeDir, fdirs, fplabs, chrs[ii])

        # get the mutation rate predicted by the model 
        mutpredict(rmutmod, xdt)

        # simulate mutations
        xdt[, "nmut" := rpois(nrow(xdt), mutRate)]
        mafdt[[ii]] <- xdt[nmut > 0L, .SD, .SDcols = c("start", "ref", "mut", "nmut")]

        rm(xdt, targetdt); gc()

    }
    nr <- sapply(mafdt, nrow)
    mafdt <- data.table::rbindlist(mafdt)

    mafdt[, "seqname" := rep(chrs, nr)]
    mafdt <- mafdt[rep(1:nrow(mafdt), nmut)]
    mafdt[, "nmut" := NULL]
    return(mafdt)

}

mutFrame <- function(rangedt, k, genomeDir, fdirs, fplabs, .chr) {

    nflank <- (k - 1) / 2
    genome <- setNames(Biostrings::readDNAStringSet(paste0(genomeDir, .chr, ".fasta")), .chr)
    sitedt <- ranges2kmerdt(rangedt$start, rangedt$end, .chr, nflank, genome)
    sitedt <- sitedt[!grepl("N", kmer)]
    pyriOrient(sitedt, isPuri(sitedt$kmer, nflank), "kmer", "kmer")
    rm(genome)

    addFeatures(fdirs, .chr, sitedt)
    xdt <- expandMuts(sitedt, nflank)
    formatFeatures(xdt, fplabs)

    return(xdt)

}

#' @export
kGet <- function(x) {

    UseMethod("kGet")

}


# these are useful accesses for glmmTMB objects
# this gets the fixed effects formula
#formula(model, fixed.only = TRUE)
# this gets strings that you can parse to make the formulas to ge the design matrices for random effects
#names(model$modelInfo$reStruc$condReStruc)