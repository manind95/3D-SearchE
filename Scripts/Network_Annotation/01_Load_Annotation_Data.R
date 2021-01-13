# Load fragments with node IDs generated from propagation outputs
DNase_Frags <- read.table("/home/maninder/Documents/PhD_Project/mESC/Annotations/Network_Frags/041220_TSS_WG_mESC_DNase_Frags.txt", sep = "\t")

# Load in annotations
Genes <- read.table("~/Documents/PhD_Project/mESC/Annotations/Genes/Coordinates.txt", sep = ",", header = T)
ChromState <- read.table("~/Documents/PhD_Project/mESC/Annotations/ChromMarks/Bed_Files/mESC_Chromatin_State.bed", sep = "\t", header = F)
  CS_posprob_files <- list.files("/home/maninder/Documents/PhD_Project/mESC/Annotations/ChromMarks/ChromHMM_mESC_Cell_Reports/POSTERIOR", ".txt", full.names = T)
P300 <- read.table("~/Documents/PhD_Project/mESC/Annotations/ChromMarks/Bed_Files/P300.bed", sep = "\t", header = F)
RNAPs2p <- read.table("~/Documents/PhD_Project/mESC/Annotations/ChromMarks/Bed_Files/RNAPII_S2P.bed", sep = "\t", header = F)
RNAPs5p <- read.table("~/Documents/PhD_Project/mESC/Annotations/ChromMarks/Bed_Files/RNAPII_S5P.bed", sep = "\t", header = F)
RNAPs7p <- read.table("~/Documents/PhD_Project/mESC/Annotations/ChromMarks/Bed_Files/RNAPII_S7P.bed", sep = "\t", header = F)
starr <- read.table("/home/maninder/Documents/PhD_Project/mESC/Annotations/Starr-seq/GSE143544_STARRseq_2iL_SL_DEseq2_mm9.tsv", sep="\t", header = TRUE)
cage <- read.table("/home/maninder/Documents/PhD_Project/mESC/Annotations/FANTOM_CAGE/mouse_permissive_enhancers_phase_1_and_2.bed", sep = "\t", header = FALSE)