##### Generate the inputs for propagation #####

# Format for props1 net.txt input 
netProp <- subset(Frags, Gene != "NG")
netProp <- netProp[,4:5]
netProp[,2] <- as.character(netProp[,2])
netProp$ID	= paste0('frag', netProp$ID)
write.table(netProp, "~/Documents/PhD_Project/mESC/Propagations/Prop1_TSS/DNase_Prop1_mESC_TSS_interactions.txt", sep = " ", quote = F, row.names = F, col.names = F)


# Load the expression matrix. 
exp_mat <- read.csv("~/Documents/PhD_Project/mESC/Propagations/Data_for_props/counttable_es_norm.csv")

# Extract lif non treated read counts and take the mean value
exp_mat <- exp_mat[,c(1,456:705)]
exp_mat$expression <- rowMeans(exp_mat[,2:251])
exp_mat <- exp_mat[,c(1,252)]

# Format for props1 bin_matrix.txt input
bin_matrix <- t(exp_mat)
write.table(bin_matrix, "~/Documents/PhD_Project/mESC/Propagations/Prop1_TSS/mESC_bin_matrix_Prop1.txt", sep = ",", quote = F, col.names = F)

# Format for props2 net.txt input

# PCHiC add chr for correct levels
#Network$baitChr <- paste("chr",Network$baitChr,sep = "")
#Network$oeChr <- paste("chr",Network$oeChr,sep = "")

# Regenerate the ChIA-PET network with the IDs and genes for each of the fragments
Network1 = Network
Frags1 = unique(Network_Frags[,1:4])

colnames(Network1) <- c("Chr", "Start", "End", "ChrOE", "StartOE", "EndOE")
test <- unique(merge(Frags1, Network1, by = c("Chr", "Start", "End"), all.y = TRUE))
colnames(Frags1) <- c("ChrOE", "StartOE", "EndOE", "IDOE")
test1 <- unique(merge(Frags1, test, by = c("ChrOE", "StartOE", "EndOE"), all.y = TRUE))
test1 <- na.omit(test1)
#HiChIP_Pilot <- test1[,c(6,7,8,9,10,3,1,2,4,5)]
netProp2 <- test1[,c(8,4)]
netProp2$ID	= paste0('frag', netProp2$ID)
netProp2$IDOE	= paste0('frag', netProp2$IDOE)
write.table(netProp2, "~/Documents/PhD_Project/mESC/Propagations/Prop2_WG/mESC_DNase_Net.txt", sep = " ", quote = F, row.names = F, col.names = F)



