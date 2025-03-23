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

pmutTargetKmers <- function(targetdt, .chr, nflank, genome) {

    xdt <- Rmutmod:::ranges2kmerdt(targetdt$start, targetdt$start, .chr, nflank, genome)

    xdt <- xdt[grepl("N", kmer) == FALSE] # remove invalid nucleotides
    pyriOrient(xdt, isPuri(xdt$kmer, nflank), "kmer", "kmer") # "mut" is assumed to already be pyrimidine oriented so we don't oriented again
    xdt[, "mut" := targetdt$mut]

    return(xdt)

}

siteTarget2pmutKmers <- function(targetdt, .chr, nflank, genome) {

    tkmerdt <- ranges2kmerdt(targetdt$start, targetdt$end, .chr, nflank, genome)
    pyriOrient(tkmerdt, isPuri(tkmerdt$kmer, nflank), "kmer", "kmer")
    tkmerdt <- tkmerdt[grepl("N", kmer) == FALSE] # remove invalid nucleotides

    # expand possible mutations
    xdt <- Rmutmod:::expandMuts(tkmerdt, nflank)

    return(xdt)

}

target2xdt <- function(targetdb, .chr, nflank, genome) {

    # identify if target is given already as possible mutations or as just sites
    tcols <- names(arrow::schema(targetdb))
    getCols <- intersect(c("start", "end", "mut"), tcols)

    # load target
    targetdt <- targetdb %>% 
        dplyr::filter(seqname == .chr) %>%
        dplyr::select(dplyr::all_of(getCols)) %>%
        dplyr::collect()

    # get kmers of each possible mutation
    xdt <- switch(
        ifelse(any(grepl("mut", getCols)), "pmuts", "sites"),
        pmuts = pmutTargetKmers(targetdt, .chr, nflank, genome),
        sites = siteTarget2pmutKmers(targetdt, .chr, nflank, genome)
    )

    return(xdt)

}

# This function takes a data frame of sequencing target sites (required columns: start, kmer), and
# a data frame of coding regions of transcripts with the right columns
# to serve as input for gtfTools::pvarCDSannotate(). The type of each possible mutation in the
# target according to their consequence in the proteins encoded by the transcripts is determined
# and all mutations other than synonymous mutations are filtered out of the target
#' @export
synonymifyTarget <- function(sitedt, cdsgtfdt, .genome) {

    # get consequence of all possible mutations in the gtf
    vardt <- gtfTools::pvarCDSannotate(cdsgtfdt, .genome)
    vardt <- vardt[type == "syn"] # conserve mutations that are synonymous in at least 1 transcript

    # get possible mutations in the target
    sitedt <- sitedt[grepl("N", kmer) == FALSE]
    sitemutdt <- Rmutmod:::expandMuts(sitedt, 0L)

    # merge possible mutations with consequence annotations and get synonymous mutations
    sitemutdt[vardt, "type" := i.type, on = setNames(c("position", "pyrimidineMut"), c("start", "mut"))]
    sitemutdt <- sitemutdt[!is.na(type), .SD, .SDcols = c("start", "mut")]

    return(sitemutdt)

}
