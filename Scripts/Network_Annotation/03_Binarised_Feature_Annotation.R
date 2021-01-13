library(GenomicRanges)
library(reshape2)


# Functions to annotate the 3C fragments and count the number of annotations per fragment
Frags_Anno <- function(Network, Annotation) {
  test <- makeGRangesFromDataFrame(Network, keep.extra.columns = T)
  test2 <- makeGRangesFromDataFrame(Annotation, keep.extra.columns = T)
  findOverlaps(test, test2) -> hits
  match_hit <- data.frame(test[queryHits(hits)] , data.frame(mcols(test2[subjectHits(hits)] )), stringsAsFactors=T)
  
  df = as(match_hit, "data.frame")
  return(df)
  
}

Count_Anno <- function(Anno, Network) {
  test <- as.data.frame(table(Anno[6]))
  colnames(test) <- c("ID", "frequency")
  test <- merge(Network, test, by="ID", all = T)
  test[is.na(test)] <- 0
  return(test)
}

# Annotate 3C frags
Frags_Anno(DNase_Frags, P300) %>%
  Count_Anno(DNase_Frags) -> test
colnames(test) <- c("ID", "Chr", "Start", "End", "P300")
Frags_Anno(DNase_Frags, RNAPs2p) %>%
  Count_Anno(DNase_Frags) -> test1
colnames(test1) <- c("ID", "Chr", "Start", "End", "RNAPs2p")
Frags_Anno(DNase_Frags, RNAPs5p) %>%
  Count_Anno(DNase_Frags) -> test2
colnames(test2) <- c("ID", "Chr", "Start", "End", "RNAPs5p")
Frags_Anno(DNase_Frags, RNAPs7p) %>%
  Count_Anno(DNase_Frags) -> test3
colnames(test3) <- c("ID", "Chr", "Start", "End", "RNAPs7p")
Frags_Anno(DNase_Frags, starr) %>%
  Count_Anno(DNase_Frags) -> test4
colnames(test3) <- c("ID", "Chr", "Start", "End", "starr")
Frags_Anno(DNase_Frags, cage) %>%
  Count_Anno(DNase_Frags) -> test5
colnames(test3) <- c("ID", "Chr", "Start", "End", "cage")
Frags_Anno(DNase_Frags, Genes) %>%
  Count_Anno(DNase_Frags) -> test6
colnames(test3) <- c("ID", "Chr", "Start", "End", "Genes")


test$RNAPs2p <- test1[,5]
test$RNAPs5p <- test2[,5]
test$RNAPs7p <- test3[,5]
test$starr <- test4[,5]
test$cage <- test5[,5]
test$Genes <- test6[,5]

DNase_Frags_Annotated <- test
