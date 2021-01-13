# Assign genes to fragments by either WG (Entrez) or TTS (Entrez_Start) coordinates
# Generate input files for propagation

# Load libraries
library(tidyverse)
library(data.table)

Genes <- read.csv("~/Documents/PhD_Project/mESC/Propagations/Data_for_props/Coordinates.txt")
Genes <- Genes[,c(2,3,4,1,5)]
Genes$Chromosome.Name <- paste("chr", Genes$Chromosome.Name, sep = "")
colnames(Genes) <- c("Chr", "Start", "End", "EnsemblID", "MGIsymbol")

# Duplicate the start coordinates to get TSS overlaps
Genes$End <- Genes$Start

# Generate a GRanges object of the unique fragments and the genes
Frags1 <- makeGRangesFromDataFrame(Network_Frags, start.field = "Start", end.field = "End", keep.extra.columns = T)
Genes <- makeGRangesFromDataFrame(unique(Genes), seqnames.field = "Chr", keep.extra.columns = T)

# Find all genic frags and append gene names
hits <- findOverlaps(Frags1, Genes) # Finds all overlaps
match_hit <- data.frame(Frags1[queryHits(hits)] , data.frame(Genes[subjectHits(hits)] ),stringsAsFactors=T)
Genic <- match_hit[,c(1,2,3,6,12)]
colnames(Genic) <- c("Chr", "Start", "End", "ID", "Gene")

# Find all non-genic frags and assign NA
non_hits <- subsetByOverlaps(Frags1, Genes, invert=TRUE) # Finds all non overlapping frags
Non_Genic <- data.frame(non_hits)
Non_Genic$Gene <- "NG"
Non_Genic <- Non_Genic[,c(1,2,3,6,7)]
colnames(Non_Genic) <- c("Chr", "Start", "End", "ID", "Gene")

# List of unique frags with id numbers and their genes
Frags <- unique(rbind(Genic, Non_Genic))
colnames(Frags) <- c("Chr", "Start", "End", "ID", "Gene")
