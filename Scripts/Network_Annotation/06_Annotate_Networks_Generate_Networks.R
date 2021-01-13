#Generate igraph objects
library(igraph)

# Load and format the network dataset
mESC_DNase_Net <- read.table("~/Documents/PhD_Project/mESC/Networks/mESC_DNaseI.txt", header = T)
mESC_DNase_Net <- subset(mESC_DNase_Net, CapC_serum_merged >= 5)
mESC_DNase_Net <- mESC_DNase_Net[,c(1,2,3,4,5,6,7)]
colnames(mESC_DNase_Net) <- c("baitChr", "baitStart", "baitEnd", "oeChr", "oeStart", "oeEnd", "Interaction_Confidence")
test <- mESC_DNase_Net[!grepl("chr1_random", mESC_DNase_Net$oeChr),]
test <- test[!grepl("chr13_random", test$oeChr),]
test <- test[!grepl("chr4_random", test$oeChr),]
test <- test[!grepl("chr8_random", test$oeChr),]
test <- test[!grepl("chr9_random", test$oeChr),]
test <- test[!grepl("chrM", test$oeChr),]
test <- test[!grepl("chrUn_random", test$oeChr),]
test <- test[!grepl("chrX_random", test$oeChr),]
mESC_DNase_Net <- test[!grepl("chrY_random", test$oeChr),]
mESC_DNase_Net$oeChr <- droplevels(mESC_DNase_Net$oeChr)

# Load and format the annotated fragment data
#mESC_DNase <- read.table("~/Documents/PhD_Project/mESC/Annotations/Network_Frags/", sep = "\t", header = T) # TODO insert the filepath and change var name

#Functions and annotated networks
mESC_MakeNets <- function(Network, Annotated_Nodes) {
  
  test <- makeGRangesFromDataFrame(Network, seqnames.field = "baitChr", start.field = "baitStart", end.field = "baitEnd", keep.extra.columns = T)
  test1 <- makeGRangesFromDataFrame(Annotated_Nodes, seqnames.field = "Chr", start.field = "Start", end.field = "End", keep.extra.columns = T)
  hits <- findOverlaps(test, test1)
  match_hit <- data.frame(test[queryHits(hits)] , data.frame(test1[subjectHits(hits)] ),stringsAsFactors=T)
  df = as(match_hit, "data.frame")
  df <- df[,c(1,2,3,6,7,8,9,15:68)]
  test <- names(df)
  test <- paste("Bait_", test, sep = "")
  colnames(df) <- test
  
  test2 <- makeGRangesFromDataFrame(df, seqnames.field = "oeChr", start.field = "oeStart", end.field = "oeEnd", keep.extra.columns = T)
  hits <- findOverlaps(test2, test1)
  match_hit <- data.frame(test2[queryHits(hits)] , data.frame(test1[subjectHits(hits)] ),stringsAsFactors=T)
  df2 = as(match_hit, "data.frame")
  df2 <- df2[,c(10,69,6,7,8,1,2,3,9,11:63,70:122)]
  test <- names(df2)
  #test <- paste("OE_", test, sep = "")
  colnames(df2) <- test 
  return(df2)
}

mESC_DNase_net_annotated <- mESC_MakeNets(mESC_DNase_Net, mESC_DNase_Frags_Annotated)
write.table(mESC_DNase_net_annotated, "/home/maninder/Documents/PhD_Project/mESC/Annotations/Network_Frags/041220_TSS_WG_mESC_DNase_Network_Annotated", sep = "\t", quote = FALSE)

# Generate igraph objects
mESC_DNase_edge_list <- mESC_DNase_net_annotated[,c(1,2,9)]

mESC_DNase_igraph <- graph.data.frame(mESC_DNase_edge_list, directed = F, vertices = mESC_DNase_Frags_Annotated) 

write.graph(mESC_DNase_igraph, "igraph_objects/mESC_DNase_igraph.graphml", format = "graphml")
















