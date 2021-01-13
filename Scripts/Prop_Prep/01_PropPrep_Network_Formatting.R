# Generates a list of unique fragments from the network and assigns each one a unique ID number

# Load libraries
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(dplyr)

# Load data
Network <- read.table("~/Documents/PhD_Project/mESC/Networks/mESC_DNaseI.txt", sep = "\t", header = T)
#Network <- read.table("~/Documents/PhD_Project/mESC/Networks/PCHiC_interaction_map.txt", sep = "\t", header = T)

# Filter for significant interactions
Network <- subset(Network, CapC_serum_merged >= 5)
#Network <- subset(Network, mESC_wt >= 5)

Network <- Network[,c(1,2,3,4,5,6)]
#Network <- Network[,c(1,2,3,6,7,8)]

# Split the interaction df and generate a single list of all the unique fragments
Network_Baits <- Network[, c(1,2,3)]
Network_OE <- Network[, c(4,5,6)]
colnames(Network_Baits) <- c("Chr", "Start", "End")
colnames(Network_OE) <- c("Chr", "Start", "End")
Network_Frags <- unique(rbind.data.frame(Network_Baits, Network_OE))

# Remove uneeded chroms
Network_Frags$Chr <- as.character(Network_Frags$Chr)
for (i in c("chrX","chrY","chr1_random","chr13_random","chr4_random","chr8_random","chr9_random","chrM","chrUn_random","chrX_random","chrY_random")) {
Network_Frags <- dplyr::filter(Network_Frags, !grepl(i, Chr))
}

# Assign unique ID's to each fragment
Network_Frags %>%
  mutate(ID = group_indices(., Chr, Start, End)) -> Network_Frags

# PCHiC add chr for correct levels
#Network_Frags$Chr <- paste("chr",Network_Frags$Chr,sep = "")

# Write frags to file
write.table(Network_Frags, file = "~/Documents/PhD_Project/mESC/Annotations/Network_Frags/TSS_WG_mESC_DNase_Frags.txt", sep = "\t", quote = F)
