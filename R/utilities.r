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

pmutTargetKmers <- function(targetdt, .chr, nflank, genome) {

    xdt <- ranges2kmerdt(targetdt$start, targetdt$start, .chr, nflank, genome)

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
    xdt <- expandMuts(tkmerdt, nflank)

    return(xdt)

}

target2xdt <- function(targetdb, .chr, nflank, genome) {

    # identify if target is given already as possible mutations or as just sites
    tcols <- names(arrow::schema(targetdb))
    getCols <- intersect(c("start", "end", "mut"), tcols)

    # load target
    seqcol <- tcols[grep("seqname|seqnames", tcols)]
    targetdt <- targetdb %>% 
        dplyr::filter(!!as.symbol(seqcol) == .chr) %>%
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
        dbcols <- names(arrow::schema(featuredb))
        seqcol <- dbcols[grep("seqname|seqnames", dbcols)]
        featuredt <- featuredb %>% 
            dplyr::filter(!!as.symbol(seqcol) == .chr) %>%
            dplyr::collect()

        lapply(dtlist, addFeature, featuredt, names(fdirs)[jj])

        jj <- jj + 1
        rm(featuredt, featuredb)

    }

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
    sitemutdt <- expandMuts(sitedt, 0L)

    # merge possible mutations with consequence annotations and get synonymous mutations
    sitemutdt[vardt, "type" := i.type, on = setNames(c("position", "pyrimidineMut"), c("start", "mut"))]
    sitemutdt <- sitemutdt[!is.na(type), .SD, .SDcols = c("start", "mut")]

    return(sitemutdt)

}

#' @export
varTargetFilter <- function(xdt, modvars, minTarget, isAggregated = FALSE) {

    badIdxs <- c()
    for (v in modvars) {

        countdt <- xdt[, list("target" = .N), by = v]
        if (isAggregated) countdt <- xdt[, list("target" = sum(nchance)), by = v]
        #print(countdt)
        badLevels <- countdt[target < minTarget][[v]]
        badIdxs <- append(badIdxs, which(xdt[[v]] %in% badLevels))

    }

    if (length(badIdxs) > 0) return(xdt[-unique(badIdxs)])
    return(xdt)

}

#' @export
dtchar2factor <- function(.dt, flevels) {

    for (jj in 1:length(flevels)) {

        v <- names(flevels)[jj]
        if (length(flevels[[jj]]) > 0) {
            .dt[, (v) := factor(as.character(get(v)), flevels[[jj]])]
        }

    }

    return()

}
