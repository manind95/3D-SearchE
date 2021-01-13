library(dplyr)
library(data.table)

# Format standardisation
Genes <- Genes[,c(2,3,4)]
Genes[,1] <- paste("chr",Genes[,1],sep = "")
colnames(Genes) <- c("Chr", "Start", "End")
ChromState <- ChromState[,c(1,2,3,4)]
colnames(ChromState) <- c("Chr", "Start", "End", "CS")
P300 <- P300[,c(1,2,3)]
colnames(P300) <- c("Chr", "Start", "End")
RNAPs2p <- RNAPs2p[,c(1,2,3)]
colnames(RNAPs2p) <- c("Chr", "Start", "End")
RNAPs5p <- RNAPs5p[,c(1,2,3)]
colnames(RNAPs5p) <- c("Chr", "Start", "End")
RNAPs7p <- RNAPs7p[,c(1,2,3)]
colnames(RNAPs7p) <- c("Chr", "Start", "End")
starr <- subset(starr, padj_SL_rep2 > 0 & padj_SL_rep1 > 0)
starr <- starr[,c(2,3,4)]
colnames(starr) <- c("Chr", "Start", "End")
cage <- cage[,1:3]
colnames(cage) <- c("chr", "start", "end")

# Filter out unwanted chroms from the fragment list
chroms <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19")

DNase_Frags$Chr <- as.character(DNase_Frags$Chr)
DNase_Frags %>% dplyr::filter(Chr %in% chroms) -> DNase_Frags

# Formatting of the posterior probabilities

test <- lapply(CS_posprob_files, fread)

chr1 <- data.frame(test[[1]])
chr2 <- data.frame(test[[12]])
chr3 <- data.frame(test[[13]])
chr4 <- data.frame(test[[14]])
chr5 <- data.frame(test[[15]])
chr6 <- data.frame(test[[16]])
chr7 <- data.frame(test[[17]])
chr8 <- data.frame(test[[18]])
chr9 <- data.frame(test[[19]])
chr10 <- data.frame(test[[2]])
chr11 <- data.frame(test[[3]])
chr12 <- data.frame(test[[4]])
chr13 <- data.frame(test[[5]])
chr14 <- data.frame(test[[6]])
chr15 <- data.frame(test[[7]])
chr16 <- data.frame(test[[8]])
chr17 <- data.frame(test[[9]])
chr18 <- data.frame(test[[10]])
chr19 <- data.frame(test[[11]])

chr1$chr <- "chr1"
chr2$chr <- "chr2"
chr3$chr <- "chr3"
chr4$chr <- "chr4"
chr5$chr <- "chr5"
chr6$chr <- "chr6"
chr7$chr <- "chr7"
chr8$chr <- "chr8"
chr9$chr <- "chr9"
chr10$chr <- "chr10"
chr11$chr <- "chr11"
chr12$chr <- "chr12"
chr13$chr <- "chr13"
chr14$chr <- "chr14"
chr15$chr <- "chr15"
chr16$chr <- "chr16"
chr17$chr <- "chr17"
chr18$chr <- "chr18"
chr19$chr <- "chr19"

chr1$start <- seq(from = 0, length.out = nrow(chr1), by=200)
chr2$start <- seq(from = 0, length.out = nrow(chr2), by=200)
chr3$start <- seq(from = 0, length.out = nrow(chr3), by=200)
chr4$start <- seq(from = 0, length.out = nrow(chr4), by=200)
chr5$start <- seq(from = 0, length.out = nrow(chr5), by=200)
chr6$start <- seq(from = 0, length.out = nrow(chr6), by=200)
chr7$start <- seq(from = 0, length.out = nrow(chr7), by=200)
chr8$start <- seq(from = 0, length.out = nrow(chr8), by=200)
chr9$start <- seq(from = 0, length.out = nrow(chr9), by=200)
chr10$start <- seq(from = 0, length.out = nrow(chr10), by=200)
chr11$start <- seq(from = 0, length.out = nrow(chr11), by=200)
chr12$start <- seq(from = 0, length.out = nrow(chr12), by=200)
chr13$start <- seq(from = 0, length.out = nrow(chr13), by=200)
chr14$start <- seq(from = 0, length.out = nrow(chr14), by=200)
chr15$start <- seq(from = 0, length.out = nrow(chr15), by=200)
chr16$start <- seq(from = 0, length.out = nrow(chr16), by=200)
chr17$start <- seq(from = 0, length.out = nrow(chr17), by=200)
chr18$start <- seq(from = 0, length.out = nrow(chr18), by=200)
chr19$start <- seq(from = 0, length.out = nrow(chr19), by=200)

chr1$end <- seq(from = 200, length.out = nrow(chr1), by=200)
chr2$end <- seq(from = 200, length.out = nrow(chr2), by=200)
chr3$end <- seq(from = 200, length.out = nrow(chr3), by=200)
chr4$end <- seq(from = 200, length.out = nrow(chr4), by=200)
chr5$end <- seq(from = 200, length.out = nrow(chr5), by=200)
chr6$end <- seq(from = 200, length.out = nrow(chr6), by=200)
chr7$end <- seq(from = 200, length.out = nrow(chr7), by=200)
chr8$end <- seq(from = 200, length.out = nrow(chr8), by=200)
chr9$end <- seq(from = 200, length.out = nrow(chr9), by=200)
chr10$end <- seq(from = 200, length.out = nrow(chr10), by=200)
chr11$end <- seq(from = 200, length.out = nrow(chr11), by=200)
chr12$end <- seq(from = 200, length.out = nrow(chr12), by=200)
chr13$end <- seq(from = 200, length.out = nrow(chr13), by=200)
chr14$end <- seq(from = 200, length.out = nrow(chr14), by=200)
chr15$end <- seq(from = 200, length.out = nrow(chr15), by=200)
chr16$end <- seq(from = 200, length.out = nrow(chr16), by=200)
chr17$end <- seq(from = 200, length.out = nrow(chr17), by=200)
chr18$end <- seq(from = 200, length.out = nrow(chr18), by=200)
chr19$end <- seq(from = 200, length.out = nrow(chr19), by=200)


probs <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19)

test <- colnames(probs[,1:20])[apply(probs[,1:20],1,which.max)]

probs$pos_prob <- test

probs %>%
  group_by(pos_prob_2 = cumsum(pos_prob != lag(pos_prob, default = dplyr::first(pos_prob)))) %>%
  mutate(E1max = max(E1)) %>%
  mutate(E2max = max(E2)) %>%
  mutate(E3max = max(E3)) %>%
  mutate(E4max = max(E4)) %>%
  mutate(E5max = max(E5)) %>%
  mutate(E6max = max(E6)) %>%
  mutate(E7max = max(E7)) %>%
  mutate(E8max = max(E8)) %>%
  mutate(E9max = max(E9)) %>%
  mutate(E10max = max(E10)) %>%
  mutate(E11max = max(E11)) %>%
  mutate(E12max = max(E12)) %>%
  mutate(E13max = max(E13)) %>%
  mutate(E14max = max(E14)) %>%
  mutate(E15max = max(E15)) %>%
  mutate(E16max = max(E16)) %>%
  mutate(E17max = max(E17)) %>%
  mutate(E18max = max(E18)) %>%
  mutate(E19max = max(E19)) %>%
  mutate(E20max = max(E20)) -> test

test <- as.data.frame(test)
test1 <- test[,c(21:24,26:45)] 

test1$prob_max <- case_when(
  test1$pos_prob == "E1" ~ test1$E1max,
  test1$pos_prob == "E2" ~ test1$E2max,
  test1$pos_prob == "E3" ~ test1$E3max,
  test1$pos_prob == "E4" ~ test1$E4max,
  test1$pos_prob == "E5" ~ test1$E5max,
  test1$pos_prob == "E6" ~ test1$E6max,
  test1$pos_prob == "E7" ~ test1$E7max,
  test1$pos_prob == "E8" ~ test1$E8max,
  test1$pos_prob == "E9" ~ test1$E9max,
  test1$pos_prob == "E10" ~ test1$E10max,
  test1$pos_prob == "E11" ~ test1$E11max,
  test1$pos_prob == "E12" ~ test1$E12max,
  test1$pos_prob == "E13" ~ test1$E13max,
  test1$pos_prob == "E14" ~ test1$E14max,
  test1$pos_prob == "E15" ~ test1$E15max,
  test1$pos_prob == "E16" ~ test1$E16max,
  test1$pos_prob == "E17" ~ test1$E17max,
  test1$pos_prob == "E18" ~ test1$E18max,
  test1$pos_prob == "E19" ~ test1$E19max,
  test1$pos_prob == "E20" ~ test1$E20max
)

CS_pos_prob <- test1

# Clean env
rm(test, test1, chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19)
