library(Rmutmod)
mafdir <- "~/projects/translateSelection/MC3/producedData/maf/"
cohort <- "KIRC"
k <- 3
targetdir <- "~/projects/translateSelection/MC3/producedData/target/"
genomedir <- "~/projects/GENCODE/release19/downloadedData/GRCh37.p13.genome/"
chrs <- c(paste0("chr", 1:22), "chrX", "chrY")
ii <- 1
.chr <- chrs[ii]
