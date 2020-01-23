#load ("~/naive_VDJ_genotyped_Dec19_Vmut3.rda") 
# H.Data is a list of data frames, where each data frame contains the final 
# genotyped data by TIgGER (after realignment using the person's own germline database as a reference)

#install.packages("dplyr")
#install.packages("stringr")
#install.packages("tidyr")
#install.packages("ggplot2")
#install.packages("ggforce")

library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggforce)

###################################

###### Determine the relative gene usage ######

full.db <- c()


for(SAMP in names(H.Data)){
  
  m = nrow(H.Data[[SAMP]])
  
  tmp.db.pre <- H.Data[[SAMP]]
  tmp.db.pre$SAMP <- SAMP
  
  # Filter out sequences with partial V-region (they usually have ambiguous gene or allele assignments)
  tmp.db.pre <- tmp.db.pre[which(tmp.db.pre$V_SEQ_LENGTH>=200),]
  
  # Remove the false positives by renaming them (to avoid skewing the gene usage by deleting them)
  tmp.db.pre$V_CALL <- gsub("\\IGHV3-9[*]01,IGHV3-9[*]01_A85C\\b", "IGHV3-9*01", tmp.db.pre$V_CALL)
  tmp.db.pre$V_CALL <- gsub("\\IGHV3-9[*]01_A85C\\b", "IGHV3-9*01", tmp.db.pre$V_CALL)
  
  tmp.db.pre$V_CALL <- gsub("\\IGHV3-9[*]01,IGHV3-9[*]01_A152G\\b", "IGHV3-9*01", tmp.db.pre$V_CALL)
  tmp.db.pre$V_CALL <- gsub("\\IGHV3-9[*]01_A152G\\b", "IGHV3-9*01", tmp.db.pre$V_CALL)
  
  tmp.db.pre$V_CALL <- gsub("\\IGHV3-11[*]01,IGHV3-11[*]01_A152G\\b", "IGHV3-11*01", tmp.db.pre$V_CALL)
  tmp.db.pre$V_CALL<- gsub("\\IGHV3-11[*]01_A152G\\b", "IGHV3-11*01", tmp.db.pre$V_CALL)
  
  tmp.db.pre$V_CALL <- gsub("\\IGHV3-11[*]05,IGHV3-11[*]05_A85C\\b", "IGHV3-11*05", tmp.db.pre$V_CALL)
  tmp.db.pre$VV_CALL <- gsub("\\IGHV3-11[*]05_A85C\\b", "IGHV3-11*05", tmp.db.pre$V_CALL)
  
  tmp.db.pre$V_CALL <- gsub("\\IGHV3-21[*]01,IGHV3-21[*]01_A152G\\b", "IGHV3-21*01", tmp.db.pre$V_CALL)
  tmp.db.pre$V_CALL <- gsub("\\IGHV3-21[*]01_A152G\\b", "IGHV3-21*01", tmp.db.pre$V_CALL)
  
  tmp.db.pre$V_CALL <- gsub("\\IGHV3-21[*]01,IGHV3-21[*]01_A85C\\b", "IGHV3-21*01", tmp.db.pre$V_CALL)
  tmp.db.pre$V_CALL <- gsub("\\IGHV3-21[*]01_A85C\\b", "IGHV3-21*01", tmp.db.pre$V_CALL)
  
  tmp.db.pre$V_CALL <- gsub("\\IGHV3-23[*]01,IGHV3-23[*]01_A85C\\b", "IGHV3-23*01", tmp.db.pre$V_CALL)
  tmp.db.pre$V_CALL <- gsub("\\IGHV3-23[*]01_A85C\\b", "IGHV3-23*01", tmp.db.pre$V_CALL)
  
  tmp.db.pre$V_CALL <- gsub("\\IGHV3-23[*]01_A152G\\b", "IGHV3-23*01", tmp.db.pre$V_CALL)
  
  tmp.db.pre$V_CALL <- gsub("\\IGHV3-30-3[*]01,IGHV3-30-3[*]01_A85C\\b", "IGHV3-30-3*01", tmp.db.pre$V_CALL)
  tmp.db.pre$V_CALL <- gsub("\\IGHV3-30-3[*]01_A85C\\b", "IGHV3-30-3*01", tmp.db.pre$V_CALL)  
  
  tmp.db.pre$V_CALL <- gsub("\\IGHV3-30-3[*]01_A152G\\b", "IGHV3-30-3*01", tmp.db.pre$V_CALL)  
  
  tmp.db.pre$V_CALL <- gsub("\\IGHV3-30[*]18,IGHV3-30[*]18_A85C\\b", "IGHV3-30*18", tmp.db.pre$V_CALL)
  tmp.db.pre$V_CALL <- gsub("\\IGHV3-30[*]18_A85C\\b", "IGHV3-30*18", tmp.db.pre$V_CALL) 
  
  tmp.db.pre$V_CALL <- gsub("\\IGHV3-33[*]01,IGHV3-33[*]01_A152G\\b", "IGHV3-33*01", tmp.db.pre$V_CALL)
  tmp.db.pre$V_CALL <- gsub("\\IGHV3-33[*]01_A152G\\b", "IGHV3-33*01", tmp.db.pre$V_CALL)
  
  tmp.db.pre$V_CALL <- gsub("\\IGHV3-33[*]01_T154G\\b", "IGHV3-33*01", tmp.db.pre$V_CALL)
  
  tmp.db.pre$V_CALL <- gsub("\\IGHV3-48[*]01,IGHV3-48[*]01_A85C\\b", "IGHV3-48*01", tmp.db.pre$V_CALL)
  tmp.db.pre$V_CALL <- gsub("\\IGHV3-48[*]01_A85C\\b", "IGHV3-48*01", tmp.db.pre$V_CALL)
  
  tmp.db.pre$V_CALL <- gsub("\\IGHV3-64D[*]06_A85C\\b", "IGHV3-64D*06", tmp.db.pre$V_CALL)
  
  # Rename selected duplicate genes and identical alleles
  tmp.db.pre$V_CALL_SINGLE <-  tmp.db.pre$V_CALL
  tmp.db.pre$V_CALL_SINGLE <- gsub("\\IGHV3-23[*]01,IGHV3-23D[*]01\\b", "IGHV3-23*01D", tmp.db.pre$V_CALL_SINGLE)
  tmp.db.pre$V_CALL_SINGLE <- gsub("\\IGHV1-69[*]01,IGHV1-69D[*]01\\b", "IGHV1-69*01D", tmp.db.pre$V_CALL_SINGLE)
  tmp.db.pre$V_CALL_SINGLE <- gsub("\\IGHV2-70[*]04,IGHV2-70D[*]04\\b", "IGHV2-70*04D", tmp.db.pre$V_CALL_SINGLE)
  tmp.db.pre$V_CALL_SINGLE <- gsub("\\IGHV3-30[*]03,IGHV3-30[*]18,IGHV3-30-5[*]01\\b", "IGHV3-30X*trip", tmp.db.pre$V_CALL_SINGLE)
  tmp.db.pre$V_CALL_SINGLE <- gsub("\\IGHV3-30[*]18,IGHV3-30-5[*]01\\b", "IGHV3-30X*doub", tmp.db.pre$V_CALL_SINGLE)
  
  # Get rid of sequences, which have ambiguous V_CALL assignment
  tmp.db.pre <- tmp.db.pre %>% rowwise %>% filter(!grepl(pattern = ',',V_CALL_SINGLE)) 
  
  # Create a list of all genes present in an individual's inferred genotype
  GENE_list <- tmp.db.pre$V_CALL_SINGLE  %>% unique() 
  GENE_list <- gsub("[*].*$","", GENE_list)
  
  # For each gene calculate the relative gene usage fraction, and keep the information about alleles in genotype
  for (gene in GENE_list){
    
    tmp.db <- tmp.db.pre  %>% filter(grepl(pattern = paste0(gene,'*'),x = V_CALL_SINGLE,fixed = T))  %>% select(V_CALL_SINGLE) %>% group_by(V_CALL_SINGLE) %>% count() %>% mutate(FRAC=n/m) 
      
    tmp.db <- separate(tmp.db, col = V_CALL_SINGLE, into = c("GENE","ALLELE"), sep = "[*]")
    
    tmp.db <- tmp.db %>% group_by(GENE) %>% mutate(ALLELE_all=paste0(unique(ALLELE),collapse=','),FRAC_SUM=sum(FRAC), N=paste0(unique(n),  collapse = ","))
    
    tmp.db$n <- NULL
    tmp.db$ALLELE <- NULL
    tmp.db$FRAC <- NULL
    
    tmp.db$SAMP <- SAMP
    
    full.db <- rbind(full.db,tmp.db)
  }
  
}



UniqueFull.db <- distinct(full.db)


###################################

###### PLOTS ######
# These plots show the relative gene usage of all individuals, 
# while also showing the different allele combinations
# in people's inferred genotypes

###### Only  IGHV1-69 ######

slctGene <- "IGHV1-69"
hv169only <- filter(UniqueFull.db, grepl(paste0("^",slctGene,"$", collapse="|"), GENE)) # ^ and $ make it find the exact pattern 
# make sure to use paste0 and not just paste, otherwise it adds a space as a separator and does not look for the exact pattern


pdf("alleleUsageFrac_HV169_20Jan20.pdf", width = 5, height = 4)
ggplot(hv169only, aes(x=ALLELE_all, y=FRAC_SUM)) +
  geom_boxplot(outlier.shape = NA, fatten=1) + 
  geom_point(position=position_dodge(width=0.4))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, size = 9, hjust = 1), axis.text.y = element_text(size = 9), strip.text = element_text(size=11), text = element_text(size=9), legend.text = element_text(size = 9)) +
  labs(title = "IGHV1-69", x="ALLELES", y="GENE USAGE FRACTION") 
dev.off()


png("alleleUsageFrac_hv169_20Jan20.png", width = 1600, height = 1200)
ggplot(hv169only, aes(x=ALLELE_all, y=FRAC_SUM)) +
  geom_boxplot(outlier.shape = NA, fatten=1) + 
  geom_point(position=position_dodge(width=0.5),aes(), size = 7)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, size = 30, hjust = 1), axis.text.y = element_text(size = 30), strip.text = element_text(size=30), text = element_text(size=30), legend.text = element_text(size = 30)) +
  labs(title = "IGHV1-69", x="ALLELES", y="GENE USAGE FRACTION") 
dev.off()



########## All remaining genes ###########


# To separate into different pages use facet_wrap_paginate from the package ggforce
# need to plot each page separately
install.packages("ggforce")
library(ggforce)

### Exclude IGHV1-69 or IGHV1-69D (allele *01 cannot be distinguished), since it has too many allele combinations and has been plotted separately using the code above
hvpattern169 <- "IGHV1-69"
hvpattern169D <- "IGHV1-69D"
UniqueNo169.db <- filter(UniqueFull.db, !grepl(paste0("^", hvpattern169,"$", collapse="|"), GENE)) # ^ and $ make it find the exact pattern # make sure to use paste0 and not just paste, otherwise it adds a space as a separator and does not look for the exact pattern
UniqueNo169.db <- filter(UniqueNo169.db, !grepl(paste0("^", hvpattern169D,"$", collapse="|"), GENE)) # optional


########## Plot all as a pdf file ##########

pdf("alleleUsageFrac_AllRemaining_boxMedian_22Jan20_4x3.pdf", width = 8, height = 10)

ggplot(UniqueNo169.db, aes(x=ALLELE_all, y=FRAC_SUM)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.5))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, size = 9, hjust = 1), axis.text.y = element_text(size = 9), strip.text = element_text(size=9)) +
  labs(x="ALLELES", y="GENE USAGE FRACTION") +
  facet_wrap_paginate(~UniqueNo169.db$GENE, nrow = 4, ncol = 3, scales = "free", page = 1)

ggplot(UniqueNo169.db, aes(x=ALLELE_all, y=FRAC_SUM)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.5))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, size = 9, hjust = 1), axis.text.y = element_text(size = 9), strip.text = element_text(size=9)) +
  labs(x="ALLELES", y="GENE USAGE FRACTION") +
  facet_wrap_paginate(~UniqueNo169.db$GENE, nrow = 4, ncol = 3, scales = "free", page = 2)

ggplot(UniqueNo169.db, aes(x=ALLELE_all, y=FRAC_SUM)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.5))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, size = 9, hjust = 1), axis.text.y = element_text(size = 9), strip.text = element_text(size=9)) +
  labs(x="ALLELES", y="GENE USAGE FRACTION") +
  facet_wrap_paginate(~UniqueNo169.db$GENE, nrow = 4, ncol = 3, scales = "free", page = 3)

ggplot(UniqueNo169.db, aes(x=ALLELE_all, y=FRAC_SUM)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.5))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, size = 9, hjust = 1), axis.text.y = element_text(size = 9), strip.text = element_text(size=9)) +
  labs(x="ALLELES", y="GENE USAGE FRACTION") +
  facet_wrap_paginate(~UniqueNo169.db$GENE, nrow = 4, ncol = 3, scales = "free", page = 4)

ggplot(UniqueNo169.db, aes(x=ALLELE_all, y=FRAC_SUM)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.5))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, size = 9, hjust = 1), axis.text.y = element_text(size = 9), strip.text = element_text(size=9)) +
  labs(x="ALLELES", y="GENE USAGE FRACTION") +
  facet_wrap_paginate(~UniqueNo169.db$GENE, nrow = 4, ncol = 3, scales = "free", page = 5)

ggplot(UniqueNo169.db, aes(x=ALLELE_all, y=FRAC_SUM)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.5))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, size = 9, hjust = 1), axis.text.y = element_text(size = 9), strip.text = element_text(size=9)) +
  labs(x="ALLELES", y="GENE USAGE FRACTION") +
  facet_wrap_paginate(~UniqueNo169.db$GENE, nrow = 4, ncol = 3, scales = "free", page = 6)

dev.off()


########## Plots as separate png files ##########

png("alleleUsageFrac_boxMedian_20Jan20_4x3_p1.png", width = 1400, height = 1800)
ggplot(UniqueNo169.db, aes(x=ALLELE_all, y=FRAC_SUM)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(size = 3)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, size = 24, hjust = 1), axis.text.y = element_text(size = 24), strip.text = element_text(size=28), text = element_text(size=26)) +
  labs(x="ALLELES", y="GENE USAGE FRACTION") +
  facet_wrap_paginate(~UniqueNo169.db$GENE, nrow = 4, ncol = 3, scales = "free", page = 1)
dev.off()

png("alleleUsageFrac_boxMedian_20Jan20_4x3_p2.png", width = 1400, height = 1800)
ggplot(UniqueNo169.db, aes(x=ALLELE_all, y=FRAC_SUM)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(size = 3)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, size = 24, hjust = 1), axis.text.y = element_text(size = 24), strip.text = element_text(size=28), text = element_text(size=26)) +
  labs(x="ALLELES", y="GENE USAGE FRACTION") +
  facet_wrap_paginate(~UniqueNo169.db$GENE, nrow = 4, ncol = 3, scales = "free", page = 2)
dev.off()

png("alleleUsageFrac_boxMedian_20Jan20_4x3_p3.png", width = 1400, height = 1800)
ggplot(UniqueNo169.db, aes(x=ALLELE_all, y=FRAC_SUM)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(size = 3)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, size = 24, hjust = 1), axis.text.y = element_text(size = 24), strip.text = element_text(size=28), text = element_text(size=26)) +
  labs(x="ALLELES", y="GENE USAGE FRACTION") +
  facet_wrap_paginate(~UniqueNo169.db$GENE, nrow = 4, ncol = 3, scales = "free", page = 3)
dev.off()

png("alleleUsageFrac_boxMedian_20Jan20_4x3_p4.png", width = 1400, height = 1800)
ggplot(UniqueNo169.db, aes(x=ALLELE_all, y=FRAC_SUM)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(size = 3)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, size = 24, hjust = 1), axis.text.y = element_text(size = 24), strip.text = element_text(size=28), text = element_text(size=26)) +
  labs(x="ALLELES", y="GENE USAGE FRACTION") +
  facet_wrap_paginate(~UniqueNo169.db$GENE, nrow = 4, ncol = 3, scales = "free", page = 4)
dev.off()

png("alleleUsageFrac_boxMedian_20Jan20_4x3_p5.png", width = 1400, height = 1800)
ggplot(UniqueNo169.db, aes(x=ALLELE_all, y=FRAC_SUM)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(size = 3)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, size = 24, hjust = 1), axis.text.y = element_text(size = 24), strip.text = element_text(size=28), text = element_text(size=26)) +
  labs(x="ALLELES", y="GENE USAGE FRACTION") +
  facet_wrap_paginate(~UniqueNo169.db$GENE, nrow = 4, ncol = 3, scales = "free", page = 5)
dev.off()

png("alleleUsageFrac_boxMedian_20Jan20_4x3_p6.png", width = 1400, height = 1800)
ggplot(UniqueNo169.db, aes(x=ALLELE_all, y=FRAC_SUM)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(size = 3)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, size = 24, hjust = 1), axis.text.y = element_text(size = 24), strip.text = element_text(size=28), text = element_text(size=26)) +
  labs(x="ALLELES", y="GENE USAGE FRACTION") +
  facet_wrap_paginate(~UniqueNo169.db$GENE, nrow = 4, ncol = 3, scales = "free", page = 6)
dev.off()

###################################


###################################
sessionInfo()
###################################
