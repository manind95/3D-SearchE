library(tidyr)
library(dplyr)
library(data.table)

# Function to find overlaps of individual chromatin states and % coverage on the fragments and join them in a df
ChromState_Anno <- function(Network, ChromState) {
  
  test <- makeGRangesFromDataFrame(Network, keep.extra.columns = T)
  test1 <- makeGRangesFromDataFrame(ChromState, keep.extra.columns = T)
  
  hits <- findOverlaps(test, test1) # Finds all overlaps
  match_hit <- data.frame(test[queryHits(hits)] , data.frame(test1[subjectHits(hits)] ),stringsAsFactors=T)
  bait_match_hit <- data.frame(test[queryHits(hits)] ,stringsAsFactors=T)
  CS_match_hit <- data.frame(test1[subjectHits(hits)]  ,stringsAsFactors=T)
  
  
  Bait_Score_CS <- makeGRangesFromDataFrame(match_hit, keep.extra.columns = T)
  x <- makeGRangesFromDataFrame(bait_match_hit, keep.extra.columns = T)
  y <- makeGRangesFromDataFrame(CS_match_hit, keep.extra.columns = T)
  overlaps <- pintersect(x,y)
  percentageOverlap <- width(overlaps) / width(x)
  match_hit <- data.frame(test[queryHits(hits)] , data.frame(mcols(test1[subjectHits(hits)] )), data.frame(percentageOverlap), stringsAsFactors=T)
  
  df = as(match_hit, "data.frame")
  return(df)
  
}
test <- ChromState_Anno(DNase_Frags, ChromState)

# Convert chromatin states to columns with percentageOverlaps as the values. I'm a motherfuckin wizard!
x=as.character(paste("seqnames+start+end+ID"))
f <- as.formula(paste(paste(x, collapse = " + "), "~ CS"))
test %>% dcast(f, value.var = "percentageOverlap", fun.aggregate = sum, fill = 0) -> test2

# Find non-overlapping nodes and add to df
test <- makeGRangesFromDataFrame(DNase_Frags, keep.extra.columns = T)
test1 <- makeGRangesFromDataFrame(ChromState, keep.extra.columns = T)
non_ov <- data.frame(test[!test %over% test1,])
non_ov <- non_ov[,c(1,2,3,6)]
non_ov$CS1 <- "LowConf"
non_ov$CS2 <- "LowConf"
non_ov$CS3 <- "LowConf"
non_ov$CS4 <- "LowConf"
non_ov$CS5 <- "LowConf"
non_ov$CS6 <- "LowConf"
non_ov$CS7 <- "LowConf"
non_ov$CS8 <- "LowConf"
non_ov$CS9 <- "LowConf"
non_ov$CS10 <- "LowConf"
non_ov$CS11 <- "LowConf"
non_ov$CS12 <- "LowConf"
non_ov$CS13 <- "LowConf"
non_ov$CS14 <- "LowConf"
non_ov$CS15 <- "LowConf"
non_ov$CS16 <- "LowConf"
non_ov$CS17 <- "LowConf"
non_ov$CS18 <- "LowConf"
non_ov$CS19 <- "LowConf"
non_ov$CS20 <- "LowConf"

colnames(test2) <- c("Chr", "Start", "End", "ID", "CS1", "CS2", "CS3", "CS4", "CS5", "CS6", "CS7", "CS8", "CS9", "CS10", "CS11", "CS12", "CS13", "CS14", "CS15", "CS16", "CS17", "CS18", "CS19", "CS20")
colnames(non_ov) <- c("Chr", "Start", "End", "ID", "CS1", "CS2", "CS3", "CS4", "CS5", "CS6", "CS7", "CS8", "CS9", "CS10", "CS11", "CS12", "CS13", "CS14", "CS15", "CS16", "CS17", "CS18", "CS19", "CS20")

DNase_CS <- rbind(non_ov, test2)