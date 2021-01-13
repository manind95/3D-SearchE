# Clean and summarise the final dataframe

#Change prop cols to numeric
mESC_DNase_Frags_Annotated$WG_n2a2_Propagation <- as.numeric(mESC_DNase_Frags_Annotated$WG_n2a2_Propagation)
mESC_DNase_Frags_Annotated$TSS_n2a2_Propagation <- as.numeric(mESC_DNase_Frags_Annotated$TSS_n2a2_Propagation)

# Most probable state per fragment
mESC_DNase_Frags_Annotated$CS_Prob <- max.col(mESC_DNase_Frags_Annotated[,c(32:51)])

# Most abundant state per fragment
mESC_DNase_Frags_Annotated$CS_Cov <- max.col(mESC_DNase_Frags_Annotated[,c(12:31)])
mESC_DNase_Frags_Annotated$CS_Cov[is.na(mESC_DNase_Frags_Annotated$CS_Cov)] <- 0

# Enhancer state coverage per fragment
mESC_DNase_Frags_Annotated$Enh_Prob <- ifelse(mESC_DNase_Frags_Annotated$E11 >= 0.95 | mESC_DNase_Frags_Annotated$E12 >= 0.95 | mESC_DNase_Frags_Annotated$E13 >= 0.95 | mESC_DNase_Frags_Annotated$E14 >= 0.95, 1, 0)

# Calculate the number of enahncer features per fragment (Collapse the RNAP data)
test <- mESC_DNase_Frags_Annotated[,c(5:10,56)]
test1 <- as.data.frame(ifelse(test >=1, 1, 0))
test2 <- rowSums(test1)
mESC_DNase_Frags_Annotated$Num_Enh_Anno <- test2

#Rename the cols
colnames(mESC_DNase_Frags_Annotated)[2:4] <- c("Chr","Start","End") 
