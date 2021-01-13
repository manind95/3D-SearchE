test1 <- makeGRangesFromDataFrame(DNase_Frags, keep.extra.columns = T)
test2 <- makeGRangesFromDataFrame(CS_pos_prob, keep.extra.columns = T)

hits <- findOverlaps(test1, test2)
match_hit <- data.frame(test1[queryHits(hits)] , data.frame(test2[subjectHits(hits)] ),stringsAsFactors=T)
df2 = as(match_hit, "data.frame")
df2 <- df2[,c(1:3,6,12:33)]
df2 <- unique(df2)

# Convert chromatin states to columns with percentageOverlaps as the values. I'm a motherfuckin wizard!
x=as.character(paste("seqnames+start+end+ID"))
f <- as.formula(paste(paste(x, collapse = " + "), "~ pos_prob"))
df2 %>% reshape2::dcast(f, value.var = "prob_max", fun.aggregate = max, fill = 0) -> test3

test3 <- test3[,c(1,2,3,4,5,16,18,19,20,21,22,23,24,6,7,8,9,10,11,12,13,14,15,17)]

DNase_CS_posprob <- test3
