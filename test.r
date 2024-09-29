# matrices
library(Rmutmod)
library(dplyr)
mafdir <- "~/projects/translateSelection/MC3/producedData/maf/"
cohort <- "PRAD"
k <- 3
targetdir <- "~/projects/translateSelection/MC3/producedData/target/"
genomedir <- "~/projects/GENCODE/release19/downloadedData/GRCh37.p13.genome/"
chrs <- paste0("chr", 1:22)
ii <- 1
.chr <- chrs[ii]
pkmers <- Rmutmod:::makePkmers(k)
fdirs <- c(
    
)
fplabs <- list(
    
)
genomePath <- "~/projects/GENCODE/release19/downloadedData/GRCh37.p13.genome/chr1.fasta"

# glm
library(Rmutmod)
library(dplyr)
mafdir <- "~/projects/translateSelection/MC3/producedData/maf/"
cohort <- "PRAD"
k <- 3
targetdir <- "~/projects/translateSelection/MC3/producedData/target/"
genomedir <- "~/projects/GENCODE/release19/downloadedData/GRCh37.p13.genome/"
chrs <- paste0("chr", 1:22)
ii <- 1
.chr <- chrs[ii]
pkmers <- Rmutmod:::makePkmers(k)
fdirs <- c(
    "f.methylation" = "~/projects/translateSelection/MC3/producedData/targetSitesFeatures/methylation/"
)
fplabs <- list(
    "f.methylation" = character(0)
)
genomePath <- "~/projects/GENCODE/release19/downloadedData/GRCh37.p13.genome/chr1.fasta"
.formula <- as.formula(
    paste0(
        "nmut ~ ",
        paste(names(fdirs), collapse = "+")
    )
)

rmutmod <- trainMutGLMs(mafdir, cohort, k, targetdir, genomedir, chrs, fdirs, fplabs, .formula)
