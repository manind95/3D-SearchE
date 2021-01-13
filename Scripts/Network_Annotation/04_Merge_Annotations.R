# Merge the annotated datsets
test <- merge(DNase_Frags_Annotated, DNase_CS, by="ID")
test <- test[,c(1:11,15:34)]
test1 <- merge(test, DNase_CS_posprob, by="ID")
test1 <- test1[,c(1:31,35:54)]
test2 <- merge(test1, DNase_Props, by="ID")
test2 <- test2[,c(1:51,55,56)]

mESC_DNase_Frags_Annotated <- test2
