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
mod2sim <- function(x, ...) {

    UseMethod("mod2sim")

}

#' @export
mod2sim.MonoMAFglmmTMB <- function(monoMAFglmmTMB, .n, pmutdt) {
    
    # load model and determine the random efects used for current mononucleotide substitution to avoid unnecessary simulations
    tmb <- modelGet(monoMAFglmmTMB)
    REs <- tmb$modelInfo$grpVar
    reList <- list()
    if (length(REs) > 0) {

        subdt <- pmutdt[, .SD, .SDcols = REs]
        reList <- lapply(subdt, unique)

    }
    
    # simulate coefficients for the linear predictor of the model
    simCoefs <- coefSim(tmb, .n, reList)

    # doing this here avoids formulas being tied to the TMB object environment
    simCoefs$cformula <- as.formula(simCoefs$cformula)
    if (length(simCoefs$ranef) > 0) {

        simCoefs$rformulas <- lapply(simCoefs$rformula, as.formula)
        simCoefs$iformulas <- lapply(simCoefs$iformula, as.formula)

    }
    
    return(simCoefs)

}


#' @export
mod2sim.MultiMAFglmmTMB <- function(multiMAFglmmTMB, .n, pmutdt) {

    tmbPaths <- modelPathsGet(multiMAFglmmTMB)
    snpTypes <- lapply(strsplit(names(tmbPaths), "[.]"), setNames, c("ref", "mut"))
    simCoefList <- list()
    for (jj in 1:length(tmbPaths)) {

        cat("\t", jj, "/", length(tmbPaths), "...\n", sep = "")

        # load model and determine the random efects used for current mononucleotide substitution to avoid unnecessary simulations
        tmb <- readRDS(tmbPaths[jj])
        REs <- tmb$modelInfo$grpVar
        reList <- list()
        if (length(REs) > 0) {

            subdt <- pmutdt[
                ref == snpTypes[[jj]]["ref"] & mut == snpTypes[[jj]]["mut"],
                .SD,
                .SDcols = REs
            ]
            reList <- lapply(subdt, unique)

        }

        # simulate coefficients for the linear predictor of the model
        simCoefList[[names(tmbPaths)[jj]]] <- coefSim(tmb, .n, reList)

        # doing this here avoids formulas being tied to the TMB object environment
        simCoefList[[names(tmbPaths)[jj]]]$cformula <- as.formula(simCoefList[[names(tmbPaths)[jj]]]$cformula)
        if (length(simCoefList[[names(tmbPaths)[jj]]]$ranef) > 0) {

            simCoefList[[names(tmbPaths)[jj]]]$rformulas <- lapply(simCoefList[[names(tmbPaths)[jj]]]$rformula, as.formula)
            simCoefList[[names(tmbPaths)[jj]]]$iformulas <- lapply(simCoefList[[names(tmbPaths)[jj]]]$iformula, as.formula)

        }
        
        rm(tmb); gc()

    }
    
    multiMAFglmmTMBsim <- new_MultiMAFglmmTMBsim(simCoefList, .n)
    return(multiMAFglmmTMBsim)

}

#' @export
coefSim <- function(tmb, .n, reList) {

    # get levels of variables
    flevels <- lapply(tmb$frame, levels)

    # if there are random effects
    .ranef <- list()
    .iformulas <- list()
    .rformulas <- list()
    if (length(reList) > 0) {
        
        # modify levels of random effect variables to be empty
        for (jj in 1:length(reList)) {

            flevels[[names(reList)[jj]]] <- character(0)

        }

        .ranef <- ranefsSim(tmb, .n, reList)
        .iformulas <- setNames(paste("~ 0 +", names(reList)), REs)
        .rformulas <- setNames(
            paste("~", sub("[|].*", "", names(tmb$modelInfo$reStruc$condReStruc))),
            names(reList)
        )

    }

    simCoefs <- new_MonoMAFglmmTMBsim(
        "fixef" = t(MASS::mvrnorm(.n, glmmTMB::fixef(tmb)$cond, vcov(tmb)$cond)),
        "ranef" = .ranef,
        "sigma" = glmmTMB::sigma(tmb),
        "cformula" = paste("~", as.character(formula(tmb, fixed.only = TRUE))[3]),
        "iformulas" = .iformulas,
        "rformulas" = .rformulas,
        "flevels" = flevels
    )

    return(simCoefs)

}

ranefsSim <- function(tmb, .n, reList) {

    # Get random effect info
    condVars <- lapply(
        glmmTMB::ranef(tmb, condVar = TRUE)$cond,
        function(.ranefs) {
            S <- attributes(.ranefs)$condVar
            S[is.na(S)] <- 0.0 # NA are zero covariances
            S 
        }
    )
    condMeans <- lapply(glmmTMB::ranef(tmb, condVar = TRUE)$cond, as.matrix, TRUE)

    # for each random effect
    reSims <- mapply(
        ranefSim,
        reList,
        condMeans,
        condVars,
        MoreArgs = list(".n" = .n),
        SIMPLIFY = FALSE
    )

    return(reSims)

}

ranefSim <- function(relevels, condMean, condVar, .n) {

    # get indices of conditional levels to simulate
    idxs <- match(relevels, rownames(condMean))

    # for each index simulate random coefficients
    L <- lapply(
        idxs,
        function(idx) t(MASS::mvrnorm(.n, condMean[idx, ], condVar[, , idx]))
    )

    return(setNames(L, rownames(condMean)[idxs]))

}


#' @export
linearPredictor <- function(x, ...) {

    UseMethod("linearPredictor", x)

}

#' @export
linearPredictor.MonoMAFglmmTMBsim <- function(simCoefs, snpdt) {

    # format fixed-effect features correctly
    flevels <- simCoefs$flevels
    for (jj in 1:length(flevels)) {

        v <- names(flevels)[jj]
        if (length(flevels[[jj]]) > 0) {
            snpdt[, (v) := factor(as.character(get(v)), flevels[[jj]])]
        }

    }

    # fixed-effects linear predictor with offset
    #FLP <- model.matrix(simCoefs$cformula, snpdt) %*% simCoefs$fixef
    FLP <- eigenMapMatMult(model.matrix(simCoefs$cformula, snpdt), simCoefs$fixef)
    mu <- FLP + snpdt$logchance

    # if there are random effects
    if (length(simCoefs$ranef) > 0) {

        # format random-effects according to those present
        reff <- names(simCoefs$ranef)
        for (jj in 1:length(reff)) {

            v <- reff[jj]
            snpdt[, (v) := factor(as.character(get(v)))]

        }

        # random-effects linear predictor
        RLP <- rePredictor(simCoefs$ranef, simCoefs$rformulas, simCoefs$iformulas, snpdt)
        mu <- mu + RLP

    }

    # add dispersion variance
    LP <- rdisp(mu, simCoefs$sigma, ncol(mu))

    # return with identifiers to later reorder matrix
    rownames(LP) <- snpdt$id
    return(LP)

}

rePredictor <- function(.ranef, rformulas, iformulas, snpdt) {

    # get levels of random effects
    rlevels <- lapply(names(.ranef), function(r) levels(snpdt[[r]]))

    # for each random effect
    Z <- list()
    for (jj in 1:length(.ranef)) {

        # get the model matrix
        M <- model.matrix(rformulas[[jj]], snpdt)

        # get the indicators matrix for the levels of the random effect
        .I <- matrix(1L, nrow = nrow(M), ncol = 1)
        if (length(rlevels[[jj]]) > 1) .I <- model.matrix(iformulas[[jj]], snpdt)

        # repeat cols of model matrix as many times as indicator column vectors
        M_blocks <- M[, rep(1:ncol(M), ncol(.I))]

        # repeat each indicator column vector as many times as columns in model matrix
        .I_blocks <- .I[, rep(1:ncol(.I), each = ncol(M))]

        # element-wise multiplication to set to 0 model matrix blocks unused in each observation
        # and setup a single matrix multiplication later
        Z[[jj]] <- M_blocks * .I_blocks
        rm(M, .I, M_blocks, .I_blocks); gc()

    }
    Z <- do.call(cbind, Z) # concatenate random effect Z matrices over the columns

    # get the random coefficient matrix for the used levels of random effects
    R <- mapply(function(L, rl) do.call(rbind, L[rl]), .ranef, rlevels, SIMPLIFY = FALSE)
    R <- do.call(rbind, R)

    # return random-effects linear predictor
    return(eigenMapMatMult(Z, R))

}

rdisp <- function(mu, .sigma, .n) {

    m <- length(mu)
    if (is.matrix(mu)) m <- nrow(mu)

    LP <- matrix(rgamma(m * .n, .sigma), nrow = m, ncol = .n) / .sigma * mu

    return(LP)

}

