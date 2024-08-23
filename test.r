library(Rmutmod)
library(dplyr)
mafdir <- "~/projects/translateSelection/MC3/producedData/maf/"
cohort <- "KIRC"
k <- 3
targetdir <- "~/projects/translateSelection/MC3/producedData/target/"
genomedir <- "~/projects/GENCODE/release19/downloadedData/GRCh37.p13.genome/"
chrs <- c(paste0("chr", 1:22), "chrX", "chrY")
ii <- 1
.chr <- chrs[ii]
pkmers <- Rmutmod:::makePkmers(k)
fdirs <- c(
    "f.tx" = "/home/mike/projects/translateSelection/MC3/producedData/targetSitesFeatures/tx"
)
fplabs <- list(
    "f.tx" = c("T", "U", "TU", "notAssignable")
)
genomePath <- "~/projects/GENCODE/release19/downloadedData/GRCh37.p13.genome/chr1.fasta"
