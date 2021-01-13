# Add on the props
prop <- fread("/home/maninder/Documents/PhD_Project/mESC/Propagations/041220_Prop_outputs/mESC_DNaseI_WG_n2a2.csv", header = FALSE, sep = ",")
prop <- transpose(prop)
prop <- as.data.frame(prop[-1,])
colnames(prop) <- c("ID", "Propagation")
prop$ID <- gsub('frag', '', prop$ID)
test <- merge(DNase_Frags, prop, by = "ID")
colnames(test) <- c("ID", "Chr","Start","End","WG_n2a2_Propagation")

# Add on the props
prop <- fread("/home/maninder/Documents/PhD_Project/mESC/Propagations/041220_Prop_outputs/mESC_DNaseI_TSS_n2a2.csv", header = FALSE, sep = ",")
prop <- transpose(prop)
prop <- as.data.frame(prop[-1,])
colnames(prop) <- c("ID", "Propagation")
prop$ID <- gsub('frag', '', prop$ID)
test1 <- merge(test, prop, by = "ID")
colnames(test1) <- c("ID", "Chr","Start","End","WG_n2a2_Propagation","TSS_n2a2_Propagation")

DNase_Props <- test1
