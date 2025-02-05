##################################
##  Landstrom Biomarker Project ##
##################################
rm(list = ls())

## Packages
library(ggplot2)
library(PCAtools)
library(pheatmap)
library(ggforce)
library(randomForest)
library(dplyr)
library(gridExtra)
library(pROC)
library(corrplot)
library(caret)
library(WriteXLS)
library(STRINGdb)
library(RColorBrewer) 
library(ggrepel)
library(glmnet)
library(caret)
library(readr)
library(readxl)
library(tidyverse)


####################################################
##  Set the folder name where the data is stored  ##
####################################################

mainfolderpath <- "XXX/ccRCC_biomarkers/" 

#############################################################################
##  Import data and check for sex and age correlation and correct for age  ##
#############################################################################

file2rd <- paste(mainfolderpath, "BiomarkerAnalysis_preAdj.csv", sep="")
Data <- read.csv(file2rd, sep = ",", head = TRUE)
rownames(Data) <- Data$SampleID
Data <- Data[,-1]
Data <- within(Data, rm("FADD")) # Remove FADD readings, see text
numcolms = ncol(Data)-14 # Separate protein list from clinical information

Data$group <- as.factor(Data$group)
Data$group2 <- as.factor(Data$group2)
Data$group_stage <- as.factor(Data$group_stage)
Data$PlateID <- as.factor(Data$PlateID)
Data$Sex <- as.factor(Data$Sex)
Data$Age <- as.numeric(Data$Age)
for (i in 1:numcolms) {Data[,i]<- as.numeric(Data[,i])}

Data_D <- Data[Data$group == "D", ]
Data_N <- Data[Data$group == "Ctrl", ]

## Linear regression with pro~age+sex in Disease (D) group
## Only 91/106 measurements are protein levels
corlist <- data.frame(matrix(NaN, numcolms, 6))
colnames(corlist) <- c("Pro", "R2","beta_age","P_age","beta_sex","P_sex")
for (i in 1:numcolms){
  corlm <- lm(Data_D[,i] ~ Data_D$Age+Data_D$Sex)
  corlist[i,1] <- colnames(Data_D)[i]
  corlist[i,2] <- summary(corlm)$r.squared
  corlist[i,3] <- coef(corlm)[[2]]
  corlist[i,4] <- summary(corlm)$coefficients[2,4]
  corlist[i,5] <- coef(corlm)[[3]]
  corlist[i,6] <- summary(corlm)$coefficients[3,4]
}
pAdj_age<- p.adjust(corlist[,4], method = "fdr")
pAdj_sex<- p.adjust(corlist[,6], method = "fdr")
corlist_D <- cbind(corlist,pAdj_age, pAdj_sex)
write.csv(corlist_D, file = paste(mainfolderpath, "cor_agesex_D.csv", sep=""))


## Linear regression with pro~age+sex in control group
corlist <- data.frame(matrix(NaN, numcolms, 6))
colnames(corlist) <- c("Pro", "R2", "beta_age", "P_age", "beta_sex", "P_sex")
for (i in 1:numcolms){
  corlm <- lm(Data_N[,i] ~ Data_N$Age + Data_N$Sex)
  corlist[i,1] <- colnames(Data_N)[i]
  corlist[i,2] <- summary(corlm)$r.squared
  corlist[i,3] <- coef(corlm)[[2]]
  corlist[i,4] <- summary(corlm)$coefficients[2,4]
  corlist[i,5] <- coef(corlm)[[3]]
  corlist[i,6] <- summary(corlm)$coefficients[3,4]
}
pAdj_age <- p.adjust(corlist[,4], method = "fdr")
pAdj_sex <- p.adjust(corlist[,6], method = "fdr")
corlist_Ctrl <- cbind(corlist,pAdj_age, pAdj_sex)
write.csv(corlist_Ctrl, file = paste(mainfolderpath, "cor_agesex_ctrl.csv", sep=""))


## Choose proteins in common
x <- corlist_D[corlist_D$pAdj_age < 0.05, ]
y <- corlist_Ctrl[corlist_Ctrl$pAdj_age < 0.05, ]
common <- merge(x, y, by.x= "Pro", by.y = "Pro")
common <- common[, -c(5,6,8,12,13,15)]

## Do age-correction
Data1 <- Data
Data1$Age <- as.numeric(Data1$Age)
Data1_D <- Data1[Data1$group == "D",]
Data1_N <- Data1[Data1$group == "Ctrl",]
for (i in 1:ncol(Data1[,-c(numcolms+1:ncol(Data1))]))
{
  if (colnames(Data1)[[i]] %in% common[, 1] == "TRUE")
  {
    n <- match(colnames(Data1)[[i]], common[, 1]) 
    beta_D <- common[n, 3]
    Data1_D[,i] <- Data1_D[,i] - beta_D * Data1_D$Age
    beta_N <- common[n,7]
    Data1_N[,i] <- Data1_N[,i] - beta_N * Data1_N$Age
    Data1 <- rbind(Data1_D, Data1_N)
    #plot(Data1_N$Age, Data1_N[,i], xlab = "age", ylab = "N_after", main = colnames(Data1_N)[[i]])
    #plot(Data_N$Age, Data_N[,i], xlab = "age", ylab = "N_before", main = colnames(Data_N)[[i]])
  }
} 
write.csv(Data1, file = paste(mainfolderpath, "BiomarkerAnalysis_ageAdj.csv", sep=""))
            

##############################################################
##  Plot the age-adjusted protein levels vs pre-adjustment  ##
##############################################################

## Load data
ageAdj_fi <- paste(mainfolderpath, "BiomarkerAnalysis_ageAdj.csv", sep="")                                      ## file path
ageAdj_dat <- read.table(ageAdj_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)            ## load to table
rownames(ageAdj_dat) <- ageAdj_dat[,1]                                                                          ## set rownames
ageAdj_dat <- ageAdj_dat[ , !(names(ageAdj_dat) %in% c("SampleID", "sampleID", ""))]                            ## delete name column

preAdj_fi <- paste(mainfolderpath, "BiomarkerAnalysis_preAdj.csv", sep="")                                      ## file path
preAdj_dat <- read.csv(preAdj_fi, sep = ",", head = TRUE)                                                       ## load to table (weird issue when using header=TRUE)
preAdj_dat <- within(preAdj_dat, rm("FADD"))                                                                    ## Remove FADD readings, see text
rownames(preAdj_dat) <- preAdj_dat$SampleID                                                                     ## set rownames
preAdj_dat <- preAdj_dat[,-1]
names(preAdj_dat) <- names(ageAdj_dat)                                                                          ## copy names from ageAdj (checked to match)

## Target proteins to do the pre/post/versus plots for
targetGenes <- common[, -c(2:9)] #("CCN4", "CD27", "CXL17", "CYR61", "EPHA2", "IGF1R", "ITGB5", "RET", "TGFBR2", "TNFRSF19", "TNFSF13", "WFDC2")
ageCorrectionPlotter <- function(geneName){
  ## before correction plot
  preAdjNPX <- preAdj_dat[[geneName]]; preAdjAGE <- preAdj_dat[["Age"]]; preAdj_frame <- data.frame(PreAdjustedNPX = preAdjNPX, Age=preAdjAGE)
  ggplot(preAdj_frame, aes(x=Age, y=PreAdjustedNPX)) + geom_point() + ggtitle(paste(geneName, " pre-correction", sep=""))
  outPath1 <- paste(mainfolderpath, "AgeCorrectionPlots/", geneName, "_preCorrection.tiff", sep="")
  ggsave(outPath1)
  
  ## post correction plot
  ageAdjNPX <- ageAdj_dat[[geneName]]; ageAdjAGE <- ageAdj_dat[["Age"]]; ageAdj_frame <- data.frame(AgeAdjustedNPX = ageAdjNPX, Age=ageAdjAGE)
  ageAdj_frame$AgeAdjustedNPX <- as.numeric(gsub(",","",ageAdj_frame$AgeAdjustedNPX))
  ggplot(ageAdj_frame, aes(x=Age, y=AgeAdjustedNPX)) + geom_point() + ggtitle(paste(geneName, " age-corrected", sep=""))
  outPath2 <- paste(mainfolderpath, "AgeCorrectionPlots/", geneName, "_ageCorrection.tiff", sep="")
  ggsave(outPath2)
  
  ## pre- vs. post-correction plot
  versus_frame <- data.frame(PreAdjustedNPX = preAdjNPX, AgeAdjustedNPX = ageAdjNPX)
  versus_frame$AgeAdjustedNPX <- as.numeric(gsub(",","",versus_frame$AgeAdjustedNPX))
  plotVersus <- ggplot(versus_frame, aes(x=PreAdjustedNPX, y=AgeAdjustedNPX)) + geom_point() + ggtitle(paste(geneName, " pre VS post-correction") )
  outPath3 <- paste(mainfolderpath, "AgeCorrectionPlots/", geneName, "_versusPlot.tiff", sep="")
  ggsave(outPath3)
}

## Run the function on all the targets
sapply(targetGenes, FUN=ageCorrectionPlotter)


###################################################################
##  Perform Wilcoxon signed-rank test and save p-values to file  ##
###################################################################

ageAdj_fi <- paste(mainfolderpath, "BiomarkerAnalysis_ageAdj.csv", sep="")                                         ## file path
ageAdj_dat <- read.table(ageAdj_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)            ## load to table
rownames(ageAdj_dat) <- ageAdj_dat[,1]                                                                          ## set rownames
ageAdj_dat <- ageAdj_dat[ , !(names(ageAdj_dat) %in% c("SampleID", "sampleID", ""))]                            ## delete name column
numcolms = ncol(ageAdj_dat)-14 # Separate protein list from clinical information

Data1_D <- ageAdj_dat[ageAdj_dat$group == "D",]
Data1_N <- ageAdj_dat[ageAdj_dat$group == "Ctrl",]

corlist <- data.frame(matrix(NaN, numcolms, 3))
colnames(corlist) <- c("Pro", "means-fold-change", "p_wilcox")
for (i in 1:numcolms){
  res <- wilcox.test(Data1_D[,i], Data1_N[,i], paired = FALSE, exact = FALSE, correct = FALSE)
  corlist[i,1] <- colnames(Data1_D)[i]
  corlist[i,2] <- mean(Data1_D[,i])/mean(Data1_N[,i])
  corlist[i,3] <- res$p.value
}
pAdj_w <- p.adjust(corlist[,3], method = "fdr")
corlist <- cbind(corlist,pAdj_w)
log2fc <- log2(corlist[,2])
corlist <- cbind(corlist,log2fc)
write.csv(corlist, file = paste(mainfolderpath, "wilcox_test_results.csv", sep=""))
            

wilcox_fi <- paste(mainfolderpath, "wilcox_test_results.csv", sep="")                                                ## file path
wilcox_dat <- read.table(wilcox_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)            ## load to table
wilcox_dat <- wilcox_dat[ , !(names(wilcox_dat) %in% c(""))]                                                    ## delete first column of numbers
wilcox_dat <- wilcox_dat[order(wilcox_dat$pAdj_w),]                                                             ## order based on p-values

## Add a column to the data frame to specify if they are UP- or DOWN- regulated
wilcox_dat$diffexpressed <- "NO"
## if log2fc > 0.2 and q-value < 0.05, set as "UP"
wilcox_dat$diffexpressed[wilcox_dat$log2fc > 0.15 & wilcox_dat$pAdj_w < 0.05] <- "UP"
## if log2fc < -0.2 and q-value < 0.05, set as "DOWN"
wilcox_dat$diffexpressed[wilcox_dat$log2fc < -0.15 & wilcox_dat$pAdj_w < 0.05] <- "DOWN"
## Add a column to label proteins on the volcano plot
wilcox_dat$label <- with(wilcox_dat, ifelse(diffexpressed=="NO", NA, Pro))
## Check
head(wilcox_dat[order(wilcox_dat$pAdj_w) & wilcox_dat$diffexpressed == 'DOWN', ])

## Create a volcano plot based on fold-change and q-values (FDR corrected p-values)
thevolcanoplot <- ggplot(data = wilcox_dat, aes(x = log2fc, y = -log10(pAdj_w),col = diffexpressed, label = label)) +
  geom_vline(xintercept = c(-0.15, 0.15), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2.5) +
  xlim(-0.75,1.5) + 
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), labels = c("Down-regulated", "Not significant", "Up-regulated")) +
  labs(color = '', x = expression("log"[2]*"(fold-change)"), y = expression("-log"[10]*"(q-value)")) +
  geom_text_repel(max.overlaps = Inf) + # To show all labels 
  theme_bw() +
  theme(axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold",size = 16),
        axis.title.y = element_text(face="bold",size = 16),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

pdf(file = paste(mainfolderpath, "wilcox_volcano.pdf", sep=""), width = 8, height = 4) # you can change the size of the output file  
thevolcanoplot # Execute the plot
dev.off() # Close the file that will contain the plot

tiff(paste(mainfolderpath, "wilcox_volcano.tiff", sep=""), units="in", width=8, height=6, res=600)
thevolcanoplot # Execute the plot
dev.off()


##################
##  Perform PCA ##
##################

## Load data
ageAdj_fi <- paste(mainfolderpath, "BiomarkerAnalysis_ageAdj.csv", sep="")                                         ## file path
ageAdj_dat <- read.table(ageAdj_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)            ## load to table
rownames(ageAdj_dat) <- ageAdj_dat[,1]                                                                          ## set rownames
ageAdj_dat <- ageAdj_dat[ , !(names(ageAdj_dat) %in% c("SampleID", "sampleID", ""))]                            ## delete name column

## split into read-data and "phenotype" data
ageAdj_values <- ageAdj_dat[1:(length(ageAdj_dat)-14)]                                                          ## remove last 14 (phenotype)
ageAdj_values <- data.frame(t(ageAdj_values))                                                                   ## flip for PCAtools
ageAdj_phenotype <- ageAdj_dat[tail(names(ageAdj_dat), 14)]                                                     ## keep only last 14 (phenotype)
rownames(ageAdj_phenotype) <- names(ageAdj_values)
ageAdj_phenotype[,'PlateID']<-factor(ageAdj_phenotype[,'PlateID'])                                              ## need integer class labels as factors

p <- pca(ageAdj_values, metadata = ageAdj_phenotype, removeVar = 0.1)                                           ## pca

tiff(paste(mainfolderpath, "PlatePCA_ageAdj.tiff", sep=""), units="in", width=5, height=5, res=600)
biplot(p, lab = NULL, colby = 'PlateID', title="Plate effect PCA", hline = 0, vline = 0, legendPosition = 'right')
dev.off()
tiff(paste(mainfolderpath, "DiseaseStatePCA_ageAdj.tiff", sep=""), units="in", width=5, height=5, res=600)
biplot(p, lab = NULL, colby = 'group', colkey = c('Ctrl' = '#00AFBB', 'D' = '#bb0c00'), title="Disease state PCA", hline = 0, vline = 0, legendPosition = 'right')
dev.off()


#####################
##  Create Heatmap ##
#####################

## Load data
ageAdj_fi <- paste(mainfolderpath, "BiomarkerAnalysis_ageAdj.csv", sep="")                                      ## file path
ageAdj_dat <- read.table(ageAdj_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)            ## load to table
rownames(ageAdj_dat) <- ageAdj_dat[,1]                                                                          ## set rownames
ageAdj_dat <- ageAdj_dat[ , !(names(ageAdj_dat) %in% c("SampleID", "sampleID", ""))]                            ## delete name column

## Metadata and Data values -- split from main frame
ageAdj_values <- ageAdj_dat[1:(length(ageAdj_dat)-14)]                                                          ## remove last 14 (phenotype)
ageAdj_values <- data.frame(t(ageAdj_values))                                                                   ## flip for PCAtools
ageAdj_phenotype <- ageAdj_dat[tail(names(ageAdj_dat), 14)]                                                     ## keep only last 14 (phenotype)
rownames(ageAdj_phenotype) <- names(ageAdj_values)
ageAdj_phenotype[,'PlateID']<-factor(ageAdj_phenotype[,'PlateID'])                                              ## need integer class labels as factors

## subset phenotype data
sub_ageAdj_phenotype <- ageAdj_phenotype[, c("Age", "Sex", "group")]

## Lists taken from wilcox test data
wilcox_fi <- paste(mainfolderpath, "wilcox_test_results.csv", sep="")                                           ## file path
wilcox_dat <- read.table(wilcox_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)            ## load to table
wilcox_dat <- wilcox_dat[ , !(names(wilcox_dat) %in% c(""))]                                                    ## delete first column of numbers
wilcox_dat <- wilcox_dat[order(wilcox_dat$pAdj_w),]                                                             ## order based on p-values

top20_proteins <- c(wilcox_dat$Pro[1:20])  ## Top20
top50_proteins <- c(wilcox_dat$Pro[1:50])  ## Top50
top70_proteins <- c(wilcox_dat$Pro[1:70])  ## Top70

top20_values <- ageAdj_values[rownames(ageAdj_values) %in% top20_proteins,]
top50_values <- ageAdj_values[rownames(ageAdj_values) %in% top50_proteins,]
top70_values <- ageAdj_values[rownames(ageAdj_values) %in% top70_proteins,]

## plot the graph(s)
plot <- pheatmap(top20_values,
                 show_rownames = T, show_colnames = T, cluster_rows = T, cluster_cols = T, clustering_distance_cols = "euclidean",
                 trace = "none", density.info = "none", annotation_col = sub_ageAdj_phenotype, annotation_names_col = T,
                 cellheight=5, cellwidth=5, fontsize=5, main='Significantly altered proteins: Top 20',
                 filename=paste(mainfolderpath, "Top20Heatmap_ageAdj.pdf", sep=""),width=20, height=20)
plot <- pheatmap(top50_values,
                 show_rownames = T, show_colnames = T, cluster_rows = T, cluster_cols = T, clustering_distance_cols = "euclidean",
                 trace = "none", density.info = "none", annotation_col = sub_ageAdj_phenotype, annotation_names_col = T,
                 cellheight=5, cellwidth=5, fontsize=5, main="Significantly altered proteins: Top 50",
                 filename=paste(mainfolderpath, "Top50Heatmap_ageAdj.pdf", sep=""),width=20, height=20)
plot <- pheatmap(top70_values,
                 show_rownames = T, show_colnames = T, cluster_rows = T, cluster_cols = T, clustering_distance_cols = "euclidean",
                 trace = "none", density.info = "none", annotation_col = sub_ageAdj_phenotype, annotation_names_col = T,
                 cellheight=5, cellwidth=5, fontsize=5, main="Significantly altered proteins: Top 70",
                 filename=paste(mainfolderpath, "Top70Heatmap_ageAdj.pdf", sep=""),width=20, height=20)


tiff(paste(mainfolderpath, "Top20Heatmap_ageAdj.tiff", sep=""), units="in", width=20, height=20, res=300)
pheatmap(top20_values, show_rownames = T, show_colnames = T, cluster_rows = T, cluster_cols = T, clustering_distance_cols = "euclidean",
         trace = "none", density.info = "none", annotation_col = sub_ageAdj_phenotype, annotation_names_col = T,
         cellheight=5, cellwidth=5, fontsize=5, main='Significantly altered proteins: Top 20')
dev.off()

tiff(paste(mainfolderpath, "Top50Heatmap_ageAdj.tiff", sep=""), units="in", width=20, height=20, res=300)
pheatmap(top50_values, show_rownames = T, show_colnames = T, cluster_rows = T, cluster_cols = T, clustering_distance_cols = "euclidean",
         trace = "none", density.info = "none", annotation_col = sub_ageAdj_phenotype, annotation_names_col = T,
         cellheight=5, cellwidth=5, fontsize=5, main="Significantly altered proteins: Top 50")
dev.off()

tiff(paste(mainfolderpath, "Top70Heatmap_ageAdj.tiff", sep=""), units="in", width=20, height=20, res=300)
pheatmap(top70_values, show_rownames = T, show_colnames = T, cluster_rows = T, cluster_cols = T, clustering_distance_cols = "euclidean",
         trace = "none", density.info = "none", annotation_col = sub_ageAdj_phenotype, annotation_names_col = T,
         cellheight=5, cellwidth=5, fontsize=5, main="Significantly altered proteins: Top 70")
dev.off()

#######################################
##  Optional: Comparative proteomics ##
#######################################

## Taking the statistically significant altered proteins and doing NPX value box plots of disease vs control
## Load data
ageAdj_fi <- paste(mainfolderpath, "BiomarkerAnalysis_ageAdj.csv", sep="")                                      ## file path
ageAdj_dat <- read.table(ageAdj_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)            ## load to table
rownames(ageAdj_dat) <- ageAdj_dat[,1]                                                                          ## set rownames
ageAdj_dat <- ageAdj_dat[ , !(names(ageAdj_dat) %in% c("SampleID", "sampleID", ""))]                            ## delete name column

## split into disease and control
disease_values <- ageAdj_dat[1:134,]
control_values <- ageAdj_dat[135:245,]

## From customer data wilcox test, top 80 significant proteins
top80_proteins <- c(wilcox_dat$Pro[1:80])

write.csv(disease_values, paste(mainfolderpath, "diseaseValues.csv", sep=""))
write.csv(control_values, paste(mainfolderpath, "controlValues.csv", sep=""))

## Using the files generated above, the boxlot information can be generated using python
## Use the provided "boxplotGenerator.py" script before running below section

## Load results from python
boxplot_fi <- paste(mainfolderpath, "boxplotData.csv", sep="")
boxplot_dat <- read.table(boxplot_fi, sep=",", header=FALSE)
names(boxplot_dat) <- c("Status","Protein","NPXValue")

## Plot boxplots
for(i in 1:4){
  p <- ggplot(data = boxplot_dat, aes(x=Status, y=NPXValue)) 
  p <- p + geom_boxplot(aes(fill=Status))
  p <- p + geom_point(aes(y=NPXValue, group=Status), position = position_dodge(width=0.75))
  p <- p + facet_wrap_paginate( ~ Protein, scales="free", nrow=5, ncol=4, page=i)
  p <- p + xlab("Status") + ylab("NPX Value") + ggtitle("Significantly altered proteins: NPX value spread")
  ggsave(paste(mainfolderpath, "BoxPlots/boxplotPage", i, ".pdf",sep=""), width=15, height=15)
}


###############################################
## Construct and test a Random Forest model  ##
###############################################

## Load data again
ageAdj_fi <- paste(mainfolderpath, "BiomarkerAnalysis_ageAdj.csv", sep="")                                      ## file path
ageAdj_dat <- read.table(ageAdj_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)            ## load to table
rownames(ageAdj_dat) <- ageAdj_dat[,1]                                                                          ## set rownames

## Lists taken from wilcox test data
wilcox_fi <- paste(mainfolderpath, "wilcox_test_results.csv", sep="")                                           ## file path
wilcox_dat <- read.table(wilcox_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)            ## load to table
wilcox_dat <- wilcox_dat[ , !(names(wilcox_dat) %in% c(""))]                                                    ## delete first column of numbers
wilcox_dat <- wilcox_dat[order(wilcox_dat$pAdj_w),]                                                             ## order based on p-values

top50_proteins <- c(wilcox_dat$Pro[1:50])
top50_proteins <- append(top50_proteins,"group")
ageAdj_top50 <- ageAdj_dat[ ,colnames(ageAdj_dat) %in% top50_proteins]
ageAdj_top50$group <- factor(ageAdj_top50$group)

## Set rng seed for reproducibility
## Make model, basic plot
set.seed(1337)
modelRF <- randomForest(group~., ageAdj_top50, ntree=100, importance=TRUE, nodesize=5, na.action=na.roughfix)
impData <- data.frame(modelRF$importance)
impData <- impData[order(impData$MeanDecreaseAccuracy, decreasing=TRUE),]
impData$Protein <- factor(rownames(impData))
rownames(impData) <- NULL

## Re-specify Protein name as relation to value (so order is retained for MeanDecreaseAccuracy)
impData$Protein = factor(impData$Protein, levels=impData[order(impData$MeanDecreaseAccuracy), "Protein"])
MDAPlot <- ggplot(impData, aes(Protein, MeanDecreaseAccuracy)) +
  geom_point(stat='identity', size=3) + 
  coord_flip() +
  labs(title="Multivariate Analysis", subtitle="Random Forest, Disease VS Control") + theme_minimal()

## Re-specify Protein name as relation to value (so order is retained for MeanDecreaseGini)
impData$Protein = factor(impData$Protein, levels=impData[order(impData$MeanDecreaseGini), "Protein"])
MDGPlot <- ggplot(impData, aes(Protein, MeanDecreaseGini)) + 
  geom_point(stat="identity", size=3) + 
  coord_flip() +
  labs(title="Multivariate Analysis", subtitle="Random Forest, Disease VS Control") + theme_minimal()

grid.arrange(MDAPlot, MDGPlot, nrow=1)
ggsave(paste(mainfolderpath, "RandomForestPlot_ageAdj_top50w.pdf", sep=""), arrangeGrob(MDAPlot, MDGPlot), width = 10, height = 9)
ggsave(paste(mainfolderpath, "RandomForestPlot_ageAdj_top50w.tiff", sep=""), arrangeGrob(MDAPlot, MDGPlot), width = 5, height = 16, dpi=600)
dev.off()

## Assess the performance of the model using cross-validation
## Split data to training and validation sets
set.seed(42)
model.data <- ageAdj_top50
model.data$class <- ifelse(model.data$group == "D", "D", "C")
model.data <- model.data[,colnames(model.data) != "group"] # Remove group column
splitSample <- createDataPartition(model.data$class, p = 0.8, list = FALSE)
training.data <- model.data[splitSample,]
validation.data <- model.data[-splitSample,]

## Run cross-validated RF
fitControl <- trainControl(
  method = "cv",
  number = 10,
  savePredictions = T,
  classProbs = T)

set.seed(13)
rfFit <- train(class ~ ., data = training.data, 
               method = "rf", 
               trControl = fitControl)

## Predictions using the best model
rfPred <- predict(rfFit, newdata = validation.data, type = "prob")
rfPred$class <- ifelse(rfPred$D > rfPred$C, "D", "C")
testPred <- validation.data$class == rfPred$class
accrcy_RF1 <- ((sum(testPred)/nrow(validation.data))*100.0)
print(accrcy_RF1)

pdf(paste(mainfolderpath, "RandomForestROC_ageAdj_top50w.pdf", sep=""))
roc <- roc(validation.data$class ~ rfPred$D, plot = TRUE, print.auc = TRUE)
dev.off()
tiff(paste(mainfolderpath, "RandomForestROC_ageAdj_top50w.tiff", sep=""), units="in", width=6, height=6, res=600)
roc <- roc(validation.data$class ~ rfPred$D, plot = TRUE, print.auc = TRUE)
dev.off()

roc.coords <- coords(roc, x = "best", ret = "all")                                                    ## Save ROC information
write.csv(roc.coords, paste(mainfolderpath, "RandomForest_ROC_metrics_ageAdj_top50w.csv", sep=""))


############################################################
## Test a Random Forest model with CV and validation sets ##
############################################################

## Load data again
ageAdj_fi <- paste(mainfolderpath, "BiomarkerAnalysis_ageAdj.csv", sep="")                                      ## file path
ageAdj_dat <- read.table(ageAdj_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)            ## load to table
rownames(ageAdj_dat) <- ageAdj_dat[,1]                                                                          ## set rownames

## Lists taken from wilcox test data
wilcox_fi <- paste(mainfolderpath, "wilcox_test_results.csv", sep="")                                           ## file path
wilcox_dat <- read.table(wilcox_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)            ## load to table
wilcox_dat <- wilcox_dat[ , !(names(wilcox_dat) %in% c(""))]                                                    ## delete first column of numbers
wilcox_dat <- wilcox_dat[order(wilcox_dat$pAdj_w),]                                                             ## order based on p-values

top50_proteins <- c(wilcox_dat$Pro[1:50])
top50_proteins <- append(top50_proteins,"group")
ageAdj_top50 <- ageAdj_dat[ ,colnames(ageAdj_dat) %in% top50_proteins]
ageAdj_top50$group <- factor(ageAdj_top50$group)

## Assess the performance of the model using cross-validation
## Split data to training and validation sets
set.seed(42)
model.data <- ageAdj_top50
model.data$class <- ifelse(model.data$group == "D", "D", "C")
model.data <- model.data[,colnames(model.data) != "group"] # Remove group column

models_RF <- data.frame(matrix(NaN, 48, 1))
accrcy_RF <- data.frame(matrix(NaN, 1, 100))
for(i in 1:100) {
    splitSample <- createDataPartition(model.data$class, p = 0.8, list = FALSE)
    training.data <- model.data[splitSample,]
    validation.data <- model.data[-splitSample,]
    
    ## Run cross-validated RF
    fitControl <- trainControl(
      method = "cv",
      number = 10,
      savePredictions = T,
      classProbs = T)
    rfFit <- train(class ~ ., data = training.data, 
                   method = "rf", 
                   trControl = fitControl)
    ## Predictions using the best model
    rfPred <- predict(rfFit, newdata = validation.data, type = "prob")
    rfPred$class <- ifelse(rfPred$D > rfPred$C, "D", "C")
    testPred <- validation.data$class == rfPred$class
    accrcy_RF[i] <- ((sum(testPred)/nrow(validation.data))*100.0)
    
    models_RF[ , paste0("trial", i)] <- rownames(validation.data)
    models_RF[ , paste0("class", i)] <- validation.data$class
    models_RF[ , paste0("pred", i)] <- rfPred$class
}

accrcy_RF <- as.matrix(accrcy_RF)
accrcy_RFm <- mean(accrcy_RF)
models_RF <- models_RF[,-1] # Remove the NaN column

write.csv(accrcy_RF,paste(mainfolderpath, "RandomForest_CV_validations_accry_top50.csv", sep=""))
write.csv(models_RF,paste(mainfolderpath, "RandomForest_CV_validations_top50.csv", sep=""))


################################################################################
## Assess the performance of the RF model with the TOP7 biomarkers 			  ##
## Do cross-validation as was done above for the top 50 significant proteins  ##
################################################################################

## Load data again
ageAdj_fi <- paste(mainfolderpath, "BiomarkerAnalysis_ageAdj.csv", sep="")                                         ## file path
ageAdj_dat <- read.table(ageAdj_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)            ## load to table
rownames(ageAdj_dat) <- ageAdj_dat[,1]                                                                          ## set rownames

set.seed(42)
model.data <- ageAdj_dat[,c("ESM1", "SDC1", "FGFBP1", "MDK", "ANXA1", "METAP2", "TFPI2","IL6", "group")]
model.data$class <- ifelse(model.data$group == "D", "D", "C")
model.data <- model.data[,colnames(model.data) != "group"] # Remove group column

models_RF_top7 <- data.frame(matrix(NaN, 48, 1))
accrcy_RF_top7 <- data.frame(matrix(NaN, 1, 100))
for(i in 1:100) {
    splitSample <- createDataPartition(model.data$class, p = 0.8, list = FALSE)
    training.data <- model.data[splitSample,]
    validation.data <- model.data[-splitSample,]
    
    # Run cross-validated RF
    fitControl <- trainControl(
      method = "cv",
      number = 10,
      savePredictions = T,
      classProbs = T)

    rfFit <- train(class ~ ., data = training.data, 
                   method = "rf", 
                   trControl = fitControl)
    
    # Predictions using the best model
    rfPred <- predict(rfFit, newdata = validation.data, type = "prob")
    rfPred$class <- ifelse(rfPred$D > rfPred$C, "D", "C")
    testPred <- validation.data$class == rfPred$class
    accrcy_RF_top7[i] <- ((sum(testPred)/nrow(validation.data))*100.0)
    
    models_RF_top7[ , paste0("trial", i)] <- rownames(validation.data)
    models_RF_top7[ , paste0("class", i)] <- validation.data$class
    models_RF_top7[ , paste0("pred", i)] <- rfPred$class
}

accrcy_RF_top7 <- as.matrix(accrcy_RF_top7)
accrcy_RF_top7m <- mean(accrcy_RF_top7)
models_RF_top7 <- models_RF_top7[,-1] # Remove the NaN column

write.csv(accrcy_RF_top7,paste(mainfolderpath, "RandomForest_CV_validations_accry_top8.csv", sep=""))
write.csv(models_RF_top7,paste(mainfolderpath, "RandomForest_CV_validations_top8.csv", sep=""))


## Do a single model for ROC curve
set.seed(42)
splitSample <- createDataPartition(model.data$class, p = 0.8, list = FALSE)
training.data <- model.data[splitSample,]
validation.data <- model.data[-splitSample,]

# Run cross-validated RF
fitControl <- trainControl(
  method = "cv",
  number = 10,
  savePredictions = T,
  classProbs = T)

rfFit <- train(class ~ ., data = training.data, 
               method = "rf", 
               trControl = fitControl)

# Predictions using the best model
rfPred <- predict(rfFit, newdata = validation.data, type = "prob")
roc <- roc(validation.data$class ~ rfPred$D, plot = TRUE, print.auc = TRUE)
pdf(paste(mainfolderpath, "RandomForestROC_top8.pdf", sep=""))
roc(validation.data$class ~ rfPred$D, plot = TRUE, print.auc = TRUE)
dev.off()

roc.coords <- coords(roc, x = "best", ret = "all")
write.csv(roc.coords,paste(mainfolderpath, "RandomForest_ROC_metrics_top8.csv", sep=""))


###########################################
##  Penalized Logistic regression models ##
###########################################

## Load data
ageAdj_fi <- paste(mainfolderpath, "BiomarkerAnalysis_ageAdj.csv", sep="")                                         ## file path
ageAdj_dat <- read.table(ageAdj_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)            ## load to table
rownames(ageAdj_dat) <- ageAdj_dat[,1]                                                                          ## set rownames
ageAdj_dat[,1] <- NULL
colnames(ageAdj_dat) <- gsub("\\.", "", colnames(ageAdj_dat))
ageAdj_dat$group <- ifelse(ageAdj_dat$group == "Ctrl", 0, 1)      ## Convert group to binary vector 0/1

## Subset to contain only columns of interest, i.e. the previously identified best markers
ageAdj_dat <- ageAdj_dat[,c("ESM1", "SDC1", "FGFBP1", "MDK", "ANXA1", "METAP2", "TFPI2", "group")]
numcols <- ncol(ageAdj_dat)

all_glm <- glmnet(as.matrix(ageAdj_dat), ageAdj_dat$group, family="binomial")
plot(all_glm, xvar='lambda', label=TRUE)

## Split data for glmmet "https://glmnet.stanford.edu/articles/glmnet.html"
set.seed(135337)
models_glm_top7 <- data.frame(matrix(NaN, 49, 1))
accrcy_glm_top7 <- data.frame(matrix(NaN, 1, 100))
for(i in 1:100) {
    splitSample <- createDataPartition(ageAdj_dat$group, p = 0.8, list = FALSE)
    training_expression <- ageAdj_dat[splitSample,]
    training_phenotype <- ageAdj_dat$group[splitSample]
    validation_expression <- ageAdj_dat[-splitSample,]
    validation_phenotype <- data.frame(Sample = rownames(validation_expression), group = validation_expression$group)
    
    ## Run cross-validated glmnet using previously defined penalisation proportion alpha
    elasticnet <- cv.glmnet(as.matrix(training_expression[,-numcols]), training_expression$group, 
                            family="binomial", nfolds=10, alpha=0.2)
    #plot(elasticnet)
    predict_validation <- predict(elasticnet, newx = as.matrix(validation_expression[,-numcols]), 
                                  s = c(elasticnet$lambda.min), type = "class")
    
    validation_group <- as.matrix(validation_phenotype[,2])
    row.names(validation_group) <- row.names(predict_validation)
    colnames(validation_group) <- colnames(predict_validation)

    testPred <- predict_validation == validation_group
    accrcy_glm_top7[i] <- (100.0-(((nrow(validation_group)-sum(testPred))/nrow(validation_group))*100.0))
    
    models_glm_top7[ , paste0("trial", i)] <- rownames(validation_expression)
    models_glm_top7[ , paste0("class", i)] <- validation_group
    models_glm_top7[ , paste0("pred", i)] <- predict_validation
}
models_glm_top7 <- models_glm_top7[,-1] # Remove the NaN column
accrcy_glm_top7 <- as.matrix(accrcy_glm_top7)
accrcy_glm_top7m <- mean(accrcy_glm_top7)

write.csv(accrcy_glm_top7,paste(mainfolderpath, "glm_CV_validations_accry_top7.csv", sep=""))
write.csv(models_glm_top7,paste(mainfolderpath, "glm_CV_validations_models_top7.csv", sep=""))

## Get regression coefficients and order by decreasing absolute coef
coefs <- coef(elasticnet)
coefs <- data.frame(coefs[-1,]) # Remove intercept
coefs <- coefs[order(abs(coefs[,1]), decreasing = T),, drop = F]

## Make ROCs by starting from the protein with largest coefficient and adding one protein at time
ord.markers <- rownames(coefs)

## Models using full data set - 7 markers
prot1Model <- glm(group ~ ANXA1, data = ageAdj_dat, family="binomial")
comb2Model <- glm(group ~ ANXA1 + METAP2, data = ageAdj_dat, family="binomial")
comb3Model <- glm(group ~ ANXA1 + METAP2 + ESM1, data = ageAdj_dat, family="binomial")
comb4Model <- glm(group ~ ANXA1 + METAP2 + ESM1 + FGFBP1, data = ageAdj_dat, family="binomial")
comb5Model <- glm(group ~ ANXA1 + METAP2 + ESM1 + FGFBP1 + MDK, data = ageAdj_dat, family="binomial")
comb6Model <- glm(group ~ ANXA1 + METAP2 + ESM1 + FGFBP1 + MDK + SDC1, data = ageAdj_dat, family="binomial")
comb7Model <- glm(group ~ ANXA1 + METAP2 + ESM1 + FGFBP1 + MDK + SDC1 + TFPI2, data = ageAdj_dat, family="binomial")

## Probabilities
prot1Prob = predict(prot1Model, newdata = ageAdj_dat, type = "response")
comb2Prob = predict(comb2Model, newdata = ageAdj_dat, type = "response")
comb3Prob = predict(comb3Model, newdata = ageAdj_dat, type = "response")
comb4Prob = predict(comb4Model, newdata = ageAdj_dat, type = "response")
comb5Prob = predict(comb5Model, newdata = ageAdj_dat, type = "response")
comb6Prob = predict(comb6Model, newdata = ageAdj_dat, type = "response")
comb7Prob = predict(comb7Model, newdata = ageAdj_dat, type = "response")

options(scipen=5)

## ROCAUC Curves
tiff(paste(mainfolderpath, "ROC_glm_top8.tiff", sep=""), units="in", width=5, height=5, res=600)
prot1ROC = roc(ageAdj_dat$group ~ prot1Prob, plot = TRUE, print.auc = FALSE, col="#e6194b",legacy.axes = TRUE)
comb2ROC = roc(ageAdj_dat$group ~ comb2Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#bcf60c',legacy.axes = TRUE)
comb3ROC = roc(ageAdj_dat$group ~ comb3Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#808000',legacy.axes = TRUE)
comb4ROC = roc(ageAdj_dat$group ~ comb4Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#3cb44b',legacy.axes = TRUE)
comb5ROC = roc(ageAdj_dat$group ~ comb5Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#fabebe',legacy.axes = TRUE)
comb6ROC = roc(ageAdj_dat$group ~ comb6Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#ffd8b1',legacy.axes = TRUE)
comb7ROC = roc(ageAdj_dat$group ~ comb7Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#469990',legacy.axes = TRUE)
dev.off()

aa <- as.data.frame(comb7Prob)
aa$class <- ifelse(aa$comb7Prob > 0.5, 1, 0)
testPred <- aa$class == ageAdj_dat$group
sum(testPred)

legend(x = 0, y = 0.5,  
       legend = c("protein1", "comb2", "comb3", "comb4", "comb5", "comb6", "comb7"),
       col = c("#e6194b", "#bcf60c", "#808000", "#3cb44b", "#fabebe", '#ffd8b1', '#469990'),
       lty = c(1, 1, 1, 1, 1, 1, 1),
       lwd = c(1, 1, 1, 1, 1, 1, 1),
       cex = 0.7)

## Prepare result table - 7 markers
res.table <- data.frame("Protein" = rownames(coefs),
                        "Coef" = coefs[,1],
                        "Combined numbers" = c("/", paste0("comb", seq(2,8))),
                        "AUC" = NA,
                        "Sen" = NA,
                        "Spe" = NA,
                        "ROC test P-value" = NA)

prob.list <- list(prot1Prob, comb2Prob, comb3Prob, comb4Prob, comb5Prob, comb6Prob, comb7Prob)
roc.list <- list(prot1ROC, comb2ROC, comb3ROC, comb4ROC, comb5ROC, comb6ROC, comb7ROC)

## AUC and other metrics for the first marker gene
AUC <- as.numeric(prot1ROC$auc)
prot1Pred <- ifelse(prot1Prob > 0.5, "1", "0")
prot1CM <- confusionMatrix(table(prot1Pred, ageAdj_dat$group))
Sen <- prot1CM$byClass["Sensitivity"]
Spe <- prot1CM$byClass["Specificity"]

res.table[1,"AUC"] <- AUC
res.table[1,"Sen"] <- Sen
res.table[1,"Spe"] <- Spe

## Add the remaining values using a loop
for(i in 2:7) {  
  res.table[i,"AUC"] <- as.numeric(roc.list[[i]]$auc)
  pred <- ifelse(prob.list[[i]] > 0.5, "1", "0")
  cm <- confusionMatrix(table(pred, ageAdj_dat$group))
  res.table[i,"Sen"] <- cm$byClass["Sensitivity"]
  res.table[i,"Spe"] <- cm$byClass["Specificity"]
}

## Add ROC test p-values
res.table[1, "ROC.test.P.value"] <- "/"
res.table[2, "ROC.test.P.value"] <- "/"
res.table[3, "ROC.test.P.value"] <- roc.test(roc.list[[2]], roc.list[[3]])$p.value
res.table[4, "ROC.test.P.value"] <- roc.test(roc.list[[3]], roc.list[[4]])$p.value
res.table[5, "ROC.test.P.value"] <- roc.test(roc.list[[4]], roc.list[[5]])$p.value
res.table[6, "ROC.test.P.value"] <- roc.test(roc.list[[5]], roc.list[[6]])$p.value
res.table[7, "ROC.test.P.value"] <- roc.test(roc.list[[6]], roc.list[[7]])$p.value

write.csv(res.table, file.path(paste(mainfolderpath, "glm_res_table_top7.csv", sep="")))

##########################################
## Models using the Top7 markers one at a time
prot1Model <- glm(group ~ ANXA1, data = ageAdj_dat, family="binomial")
prot2Model <- glm(group ~ METAP2, data = ageAdj_dat, family="binomial")
prot3Model <- glm(group ~ ESM1, data = ageAdj_dat, family="binomial")
prot4Model <- glm(group ~ FGFBP1, data = ageAdj_dat, family="binomial")
prot5Model <- glm(group ~ MDK, data = ageAdj_dat, family="binomial")
prot6Model <- glm(group ~ SDC1, data = ageAdj_dat, family="binomial")
prot7Model <- glm(group ~ TFPI2, data = ageAdj_dat, family="binomial")

## Probabilities
prot1Prob = predict(prot1Model, newdata = ageAdj_dat, type = "response")
prot2Prob = predict(prot2Model, newdata = ageAdj_dat, type = "response")
prot3Prob = predict(prot3Model, newdata = ageAdj_dat, type = "response")
prot4Prob = predict(prot4Model, newdata = ageAdj_dat, type = "response")
prot5Prob = predict(prot5Model, newdata = ageAdj_dat, type = "response")
prot6Prob = predict(prot6Model, newdata = ageAdj_dat, type = "response")
prot7Prob = predict(prot7Model, newdata = ageAdj_dat, type = "response")

## ROCAUC Curves
tiff(paste(mainfolderpath, "ROC_glm_top8_singles.tiff", sep=""), units="in", width=5, height=5, res=600)
prot1ROC = roc(ageAdj_dat$group ~ prot1Prob, plot = TRUE, print.auc = FALSE, col="#e6194b",legacy.axes = TRUE)
prot2ROC = roc(ageAdj_dat$group ~ prot2Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#bcf60c',legacy.axes = TRUE)
prot3ROC = roc(ageAdj_dat$group ~ prot3Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#808000',legacy.axes = TRUE)
prot4ROC = roc(ageAdj_dat$group ~ prot4Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#3cb44b',legacy.axes = TRUE)
prot5ROC = roc(ageAdj_dat$group ~ prot5Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#fabebe',legacy.axes = TRUE)
prot6ROC = roc(ageAdj_dat$group ~ prot6Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#ffd8b1',legacy.axes = TRUE)
prot7ROC = roc(ageAdj_dat$group ~ prot7Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#469990',legacy.axes = TRUE)
dev.off()

aa <- as.data.frame(prot1Prob)
aa$class <- ifelse(aa$prot1Prob > 0.5, 1, 0)
testPred <- aa$class == ageAdj_dat$group
sum(testPred)

legend(x = 0, y = 0.5,  
       legend = c("ANXA1", "METAP2", "ESM1", "SDC1", "FGFBP1", "MDK", "TFPI2"),
       col = c("#e6194b", "#bcf60c", "#808000", "#3cb44b", "#fabebe", '#ffd8b1', '#469990'),
       lty = c(1, 1, 1, 1, 1, 1, 1),
       lwd = c(1, 1, 1, 1, 1, 1, 1),
       cex = 0.7)

## Prepare result table - 7 markers
res.table_ind <- data.frame("Protein" = rownames(coefs),
                        "Coef" = coefs[,1],
                        "Combined numbers" = c("/", paste0("prot", seq(2,7))),
                        "AUC" = NA,
                        "Sen" = NA,
                        "Spe" = NA,
                        "ROC test P-value" = NA)

prob.list <- list(prot1Prob, prot2Prob, prot3Prob, prot4Prob, prot5Prob, prot6Prob, prot7Prob)
roc.list <- list(prot1ROC, prot2ROC, prot3ROC, prot4ROC, prot5ROC, prot6ROC, prot7ROC)

## AUC and other metrics for the first marker gene
AUC <- as.numeric(prot1ROC$auc)
prot1Pred <- ifelse(prot1Prob > 0.5, "1", "0")
prot1CM <- confusionMatrix(table(prot1Pred, ageAdj_dat$group))
Sen <- prot1CM$byClass["Sensitivity"]
Spe <- prot1CM$byClass["Specificity"]

res.table_ind[1,"AUC"] <- AUC
res.table_ind[1,"Sen"] <- Sen
res.table_ind[1,"Spe"] <- Spe

## Add the remaining values using a loop
for(i in 2:7) {  
  res.table_ind[i,"AUC"] <- as.numeric(roc.list[[i]]$auc)
  pred <- ifelse(prob.list[[i]] > 0.5, "1", "0")
  cm <- confusionMatrix(table(pred, ageAdj_dat$group))
  res.table_ind[i,"Sen"] <- cm$byClass["Sensitivity"]
  res.table_ind[i,"Spe"] <- cm$byClass["Specificity"]
}

write.csv(res.table_ind, file.path(paste(mainfolderpath, "glm_res_table_top7_singles.csv", sep="")))

##########################################################
##  Penalized Logistic regression with 7 random markers ##
##########################################################

set.seed(135337)
models_glm_rdm7_genes <- data.frame(matrix(NaN, 1, 7))
for(j in 1:10) {
  ## Load data
  ageAdj_fi <- paste(mainfolderpath, "BiomarkerAnalysis_ageAdj.csv", sep="")                                      ## file path
  ageAdj_dat <- read.table(ageAdj_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)            ## load to table
  rownames(ageAdj_dat) <- ageAdj_dat[,1]                                                                          ## set rownames
  ageAdj_dat[,1] <- NULL
  colnames(ageAdj_dat) <- gsub("\\.", "", colnames(ageAdj_dat))
  ageAdj_dat$group <- ifelse(ageAdj_dat$group == "Ctrl", 0, 1)      ## Convert group to binary vector 0/1
  
  ## Subset to contain only columns of interest, here, just random 7
  ageAdj_dat_rdm7 <- sample(ageAdj_dat[,1:91], 7)  ##ageAdj_dat <- ageAdj_dat[,c("IFNGR1", "CEACAM1", "CD207", "AREG", "FURIN", "MUC16", "CD27", "group")]
  ageAdj_dat_rdm7$group <- ageAdj_dat$group
  ageAdj_dat <- ageAdj_dat_rdm7
  
  numcols <- ncol(ageAdj_dat)
  models_glm_rdm7_genes[paste0("trial", j), ] <- colnames(ageAdj_dat_rdm7)[-8]
  
  all_glm <- glmnet(as.matrix(ageAdj_dat), ageAdj_dat$group, family="binomial")
  #plot(all_glm, xvar='lambda', label=TRUE)
  
  ## Split data for glmmet
  accrcy_glm_rdm7 <- data.frame(matrix(NaN, 1, 100))
  models_glm_rdm7 <- data.frame(matrix(NaN, 49, 1))
  for(i in 1:100) {
    splitSample <- createDataPartition(ageAdj_dat$group, p = 0.8, list = FALSE)
    training_expression <- ageAdj_dat[splitSample,]
    training_phenotype <- ageAdj_dat$group[splitSample]
    validation_expression <- ageAdj_dat[-splitSample,]
    validation_phenotype <- data.frame(Sample = rownames(validation_expression), group = validation_expression$group)
    
    ## Run cross-validated glmnet using previously defined penalisation proportion alpha
    elasticnet <- cv.glmnet(as.matrix(training_expression[,-numcols]), training_expression$group, 
                            family="binomial", nfolds=10, alpha=0.2)
    #plot(elasticnet)
    predict_validation <- predict(elasticnet, newx = as.matrix(validation_expression[,-numcols]), 
                                  s = c(elasticnet$lambda.min), type = "class")
    
    validation_group <- as.matrix(validation_phenotype[,2])
    row.names(validation_group) <- row.names(predict_validation)
    colnames(validation_group) <- colnames(predict_validation)
    
    testPred <- predict_validation == validation_group
    accrcy_glm_rdm7[i] <- (100.0-(((nrow(validation_group)-sum(testPred))/nrow(validation_group))*100.0))
    
    models_glm_rdm7[ , paste0("trial", i)] <- rownames(validation_expression)
    models_glm_rdm7[ , paste0("class", i)] <- validation_group
    models_glm_rdm7[ , paste0("pred", i)] <- predict_validation
  }
  
  accrcy_glm_rdm7 <- as.matrix(accrcy_glm_rdm7)
  accrcy_glm_rdm7m <- mean(accrcy_glm_rdm7)
  models_glm_rdm7 <- models_glm_rdm7[,-1] # Remove the NaN column
  
  write.csv(accrcy_glm_rdm7,paste(mainfolderpath, "glm_CV_validations_accry_rdm7_", j, ".csv", sep=""))
  write.csv(models_glm_rdm7,paste(mainfolderpath, "glm_CV_validations_models_rdm7_", j, ".csv", sep=""))
  
}

models_glm_rdm7_genes <- models_glm_rdm7_genes[-1,]
write.csv(models_glm_rdm7_genes,paste(mainfolderpath, "glm_CV_validations_models_rdm7_genes.csv", sep=""))

##########################################
## Get another instance of random 7, bottom 7, or non-significant 7 and plot ROC curves and stats

set.seed(141)
ageAdj_fi <- paste(mainfolderpath, "BiomarkerAnalysis_ageAdj.csv", sep="")                                         ## file path
ageAdj_dat <- read.table(ageAdj_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)            ## load to table
rownames(ageAdj_dat) <- ageAdj_dat[,1]                                                                          ## set rownames
ageAdj_dat[,1] <- NULL
colnames(ageAdj_dat) <- gsub("\\.", "", colnames(ageAdj_dat))
ageAdj_dat$group <- ifelse(ageAdj_dat$group == "Ctrl", 0, 1)      ## Convert group to binary vector 0/1

wilcox_fi <- paste(mainfolderpath, "wilcox_test_results.csv", sep="")                                                ## file path
wilcox_dat <- read.table(wilcox_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)            ## load to table
wilcox_dat <- wilcox_dat[ , !(names(wilcox_dat) %in% c(""))]                                ## delete first column of numbers
wilcox_dat <- wilcox_dat[order(wilcox_dat$pAdj_w),] ## order based on p-values

## Subset to contain only bottom 7 of top50
bot7_proteins <- c(wilcox_dat$Pro[44:50])
bot7_proteins <- append(bot7_proteins,"group")
ageAdj_bot7 <- ageAdj_dat[ ,colnames(ageAdj_dat) %in% bot7_proteins]
ageAdj_bot7$group <- factor(ageAdj_bot7$group)

## Subset to contain only bottom 7 of all
botbot7_proteins <- c(wilcox_dat$Pro[85:91])
botbot7_proteins <- append(botbot7_proteins,"group")
ageAdj_botbot7 <- ageAdj_dat[ ,colnames(ageAdj_dat) %in% botbot7_proteins]
ageAdj_botbot7$group <- factor(ageAdj_botbot7$group)

## Subset to contain only columns of interest, here, just random 7
ageAdj_dat_rdm7 <- sample(ageAdj_dat[,1:91], 7)
ageAdj_dat_rdm7$group <- ageAdj_dat$group
ageAdj_dat_rdm7$group <- factor(ageAdj_dat_rdm7$group)

#################################################################
## Choose one below and un-comment as necessary
## random 7 from all data
ageAdj_dat <- ageAdj_dat_rdm7

## or bottom 7 from top50
#ageAdj_dat <- ageAdj_bot7

## or bottom 7, non-significant ones from all
#ageAdj_dat <- ageAdj_botbot7
#################################################################

numcols <- ncol(ageAdj_dat)
models_glm_rdm7_genes_plt <- data.frame(matrix(NaN, 1, 7))
models_glm_rdm7_genes_plt["plt", ] <- colnames(ageAdj_dat)[-8]
models_glm_rdm7_genes_plt <- models_glm_rdm7_genes_plt[-1,]

set.seed(142)
## Split data for glmmet
splitSample <- createDataPartition(ageAdj_dat$group, p = 0.8, list = FALSE)
training_expression <- ageAdj_dat[splitSample,]
training_phenotype <- ageAdj_dat$group[splitSample]
validation_expression <- ageAdj_dat[-splitSample,]
validation_phenotype <- data.frame(Sample = rownames(validation_expression), group = validation_expression$group)

## Run cross-validated glmnet using previously defined penalisation proportion alpha
elasticnet <- cv.glmnet(as.matrix(training_expression[,-numcols]), training_expression$group, 
                        family="binomial", nfolds=10, alpha=0.2)
#plot(elasticnet)
predict_validation <- predict(elasticnet, newx = as.matrix(validation_expression[,-numcols]), 
                            s = c(elasticnet$lambda.min), type = "class")

validation_group <- as.matrix(validation_phenotype[,2])
row.names(validation_group) <- row.names(predict_validation)
colnames(validation_group) <- colnames(predict_validation)

testPred <- predict_validation == validation_group
accrcy_glm_rdm7 <- (100.0-(((nrow(validation_group)-sum(testPred))/nrow(validation_group))*100.0))

models_glm_rdm7 <- data.frame(matrix(NaN, 48, 1))
models_glm_rdm7[ , "trial"] <- rownames(validation_expression)
models_glm_rdm7[ , "class"] <- validation_group
models_glm_rdm7[ , "pred"] <- predict_validation
models_glm_rdm7 <- models_glm_rdm7[,-1] # Remove the NaN column

#################################################################
write.csv(models_glm_rdm7_genes_plt,paste(mainfolderpath, "glm_CV_validations_models_rdm7_genes_plt.csv", sep=""))
write.csv(accrcy_glm_rdm7,paste(mainfolderpath, "glm_CV_validations_accry_rdm7_plt.csv", sep=""))
write.csv(models_glm_rdm7,paste(mainfolderpath, "glm_CV_validations_models_rdm7_plt.csv", sep=""))

#write.csv(models_glm_rdm7_genes_plt,paste(mainfolderpath, "glm_CV_validations_models_bot7_genes_plt.csv", sep=""))     ## un-comment as necessary
#write.csv(accrcy_glm_rdm7,paste(mainfolderpath, "glm_CV_validations_accry_bot7_plt.csv", sep=""))
#write.csv(models_glm_rdm7,paste(mainfolderpath, "glm_CV_validations_models_bot7_plt.csv", sep=""))

#write.csv(models_glm_rdm7_genes_plt,paste(mainfolderpath, "glm_CV_validations_models_botbot7_genes_plt.csv", sep=""))  ## un-comment as necessary
#write.csv(accrcy_glm_rdm7,paste(mainfolderpath, "glm_CV_validations_accry_botbot7_plt.csv", sep=""))
#write.csv(models_glm_rdm7,paste(mainfolderpath, "glm_CV_validations_models_botbot7_plt.csv", sep=""))
#################################################################

## Get regression coefficients and order by decreasing absolute coef
coefs <- coef(elasticnet)
coefs <- data.frame(coefs[-1,]) # Remove intercept
coefs <- coefs[order(abs(coefs[,1]), decreasing = T),, drop = F]

## Make ROCs by starting from the protein with largest coefficient and adding one protein at time
ord.markers <- rownames(coefs)

## Models using random 7 markers
prot1Model <- glm(group ~ CXCL13, data = ageAdj_dat, family="binomial")
comb2Model <- glm(group ~ CXCL13 + VEGFA, data = ageAdj_dat, family="binomial")
comb3Model <- glm(group ~ CXCL13 + VEGFA + PODXL, data = ageAdj_dat, family="binomial")
comb4Model <- glm(group ~ CXCL13 + VEGFA + PODXL + TRAIL, data = ageAdj_dat, family="binomial")
comb5Model <- glm(group ~ CXCL13 + VEGFA + PODXL + TRAIL + SMAD5, data = ageAdj_dat, family="binomial")
comb6Model <- glm(group ~ CXCL13 + VEGFA + PODXL + TRAIL + SMAD5 + GPC1, data = ageAdj_dat, family="binomial")
comb7Model <- glm(group ~ CXCL13 + VEGFA + PODXL + TRAIL + SMAD5 + GPC1 + TGFA, data = ageAdj_dat, family="binomial")

## Models using bottom (of top 50) 7 markers
prot1Model <- glm(group ~ CEACAM1, data = ageAdj_dat, family="binomial")
comb2Model <- glm(group ~ CEACAM1 + FASLG, data = ageAdj_dat, family="binomial")
comb3Model <- glm(group ~ CEACAM1 + FASLG + CD48, data = ageAdj_dat, family="binomial")
comb4Model <- glm(group ~ CEACAM1 + FASLG + CD48 + CD160, data = ageAdj_dat, family="binomial")
comb5Model <- glm(group ~ CEACAM1 + FASLG + CD48 + CD160 + MUC16, data = ageAdj_dat, family="binomial")
comb6Model <- glm(group ~ CEACAM1 + FASLG + CD48 + CD160 + MUC16 + CDKN1A, data = ageAdj_dat, family="binomial")
comb7Model <- glm(group ~ CEACAM1 + FASLG + CD48 + CD160 + MUC16 + CDKN1A + ICOSLG, data = ageAdj_dat, family="binomial")

## Models using bottom/non-significant (of all markers, based on wilcox test) 7 markers
prot1Model <- glm(group ~ LYPD3, data = ageAdj_dat, family="binomial")
comb2Model <- glm(group ~ LYPD3 + SCF, data = ageAdj_dat, family="binomial")
comb3Model <- glm(group ~ LYPD3 + SCF + SPARC, data = ageAdj_dat, family="binomial")
comb4Model <- glm(group ~ LYPD3 + SCF + SPARC + GZMH, data = ageAdj_dat, family="binomial")
comb5Model <- glm(group ~ LYPD3 + SCF + SPARC + GZMH + EPHA2, data = ageAdj_dat, family="binomial")
comb6Model <- glm(group ~ LYPD3 + SCF + SPARC + GZMH + EPHA2 + CRNN, data = ageAdj_dat, family="binomial")
comb7Model <- glm(group ~ LYPD3 + SCF + SPARC + GZMH + EPHA2 + CRNN + CXL17, data = ageAdj_dat, family="binomial")

## Probabilities
prot1Prob = predict(prot1Model, newdata = ageAdj_dat, type = "response")
comb2Prob = predict(comb2Model, newdata = ageAdj_dat, type = "response")
comb3Prob = predict(comb3Model, newdata = ageAdj_dat, type = "response")
comb4Prob = predict(comb4Model, newdata = ageAdj_dat, type = "response")
comb5Prob = predict(comb5Model, newdata = ageAdj_dat, type = "response")
comb6Prob = predict(comb6Model, newdata = ageAdj_dat, type = "response")
comb7Prob = predict(comb7Model, newdata = ageAdj_dat, type = "response")

## ROCAUC Curves
tiff(paste(mainfolderpath, "ROC_glm_botbot7_plt.tiff", sep=""), units="in", width=5, height=5, res=600)
prot1ROC = roc(ageAdj_dat$group ~ prot1Prob, plot = TRUE, print.auc = FALSE, col="#e6194b",legacy.axes = TRUE)
comb2ROC = roc(ageAdj_dat$group ~ comb2Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#bcf60c',legacy.axes = TRUE)
comb3ROC = roc(ageAdj_dat$group ~ comb3Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#808000',legacy.axes = TRUE)
comb4ROC = roc(ageAdj_dat$group ~ comb4Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#3cb44b',legacy.axes = TRUE)
comb5ROC = roc(ageAdj_dat$group ~ comb5Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#fabebe',legacy.axes = TRUE)
comb6ROC = roc(ageAdj_dat$group ~ comb6Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#ffd8b1',legacy.axes = TRUE)
comb7ROC = roc(ageAdj_dat$group ~ comb7Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#469990',legacy.axes = TRUE)
dev.off()

aa <- as.data.frame(comb7Prob)
aa$class <- ifelse(aa$comb7Prob > 0.5, 1, 0)
testPred <- aa$class == ageAdj_dat$group
sum(testPred)

legend(x = 0, y = 0.5,  
       legend = c("protein1", "comb2", "comb3", "comb4", "comb5", "comb6", "comb7"),
       col = c("#e6194b", "#bcf60c", "#808000", "#3cb44b", "#fabebe", '#ffd8b1', '#469990'),
       lty = c(1, 1, 1, 1, 1, 1, 1),
       lwd = c(1, 1, 1, 1, 1, 1, 1),
       cex = 0.7)

## Prepare result table - random 7 markers
res.table <- data.frame("Protein" = rownames(coefs),
                        "Coef" = coefs[,1],
                        "Combined numbers" = c("/", paste0("comb", seq(2,7))),
                        "AUC" = NA,
                        "Sen" = NA,
                        "Spe" = NA,
                        "ROC test P-value" = NA)

prob.list <- list(prot1Prob, comb2Prob, comb3Prob, comb4Prob, comb5Prob, comb6Prob, comb7Prob)
roc.list <- list(prot1ROC, comb2ROC, comb3ROC, comb4ROC, comb5ROC, comb6ROC, comb7ROC)

## AUC and other metrics for the first marker gene
AUC <- as.numeric(prot1ROC$auc)
prot1Pred <- ifelse(prot1Prob > 0.5, "1", "0")
prot1CM <- confusionMatrix(table(prot1Pred, ageAdj_dat$group))
Sen <- prot1CM$byClass["Sensitivity"]
Spe <- prot1CM$byClass["Specificity"]

res.table[1,"AUC"] <- AUC
res.table[1,"Sen"] <- Sen
res.table[1,"Spe"] <- Spe

## Add the remaining values using a loop
for(i in 2:7) {  
  res.table[i,"AUC"] <- as.numeric(roc.list[[i]]$auc)
  pred <- ifelse(prob.list[[i]] > 0.5, "1", "0")
  cm <- confusionMatrix(table(pred, ageAdj_dat$group))
  res.table[i,"Sen"] <- cm$byClass["Sensitivity"]
  res.table[i,"Spe"] <- cm$byClass["Specificity"]
}

## Add ROC test p-values
res.table[1, "ROC.test.P.value"] <- "/"
res.table[2, "ROC.test.P.value"] <- "/"
res.table[3, "ROC.test.P.value"] <- roc.test(roc.list[[2]], roc.list[[3]])$p.value
res.table[4, "ROC.test.P.value"] <- roc.test(roc.list[[3]], roc.list[[4]])$p.value
res.table[5, "ROC.test.P.value"] <- roc.test(roc.list[[4]], roc.list[[5]])$p.value
res.table[6, "ROC.test.P.value"] <- roc.test(roc.list[[5]], roc.list[[6]])$p.value
res.table[7, "ROC.test.P.value"] <- roc.test(roc.list[[6]], roc.list[[7]])$p.value

#################################################################
write.csv(res.table, file.path(paste(mainfolderpath, "glm_res_table_rdm7_plt.csv", sep="")))
#write.csv(res.table, file.path(paste(mainfolderpath, "glm_res_table_bot7_plt.csv", sep="")))         ## un-comment as necessary
#write.csv(res.table, file.path(paste(mainfolderpath, "glm_res_table_botbot7_plt.csv", sep="")))      ## un-comment as necessary
#################################################################

ageAdj_dat <- ageAdj_bot7      ## un-comment as necessary
#ageAdj_dat <- ageAdj_botbot7   ## un-comment as necessary
numcols <- ncol(ageAdj_dat)

## Split data for glmmet
accrcy_glm_rdm7 <- data.frame(matrix(NaN, 1, 100))
models_glm_rdm7 <- data.frame(matrix(NaN, 48, 1))
for(i in 1:100) {
  splitSample <- createDataPartition(ageAdj_dat$group, p = 0.8, list = FALSE)
  training_expression <- ageAdj_dat[splitSample,]
  training_phenotype <- ageAdj_dat$group[splitSample]
  validation_expression <- ageAdj_dat[-splitSample,]
  validation_phenotype <- data.frame(Sample = rownames(validation_expression), group = validation_expression$group)
  
  ## Run cross-validated glmnet using previously defined penalisation proportion alpha
  elasticnet <- cv.glmnet(as.matrix(training_expression[,-numcols]), training_expression$group, 
                          family="binomial", nfolds=10, alpha=0.2)
  #plot(elasticnet)
  predict_validation <- predict(elasticnet, newx = as.matrix(validation_expression[,-numcols]), 
                                s = c(elasticnet$lambda.min), type = "class")
  
  validation_group <- as.matrix(validation_phenotype[,2])
  row.names(validation_group) <- row.names(predict_validation)
  colnames(validation_group) <- colnames(predict_validation)
  
  testPred <- predict_validation == validation_group
  accrcy_glm_rdm7[i] <- (100.0-(((nrow(validation_group)-sum(testPred))/nrow(validation_group))*100.0))
  
  models_glm_rdm7[ , paste0("trial", i)] <- rownames(validation_expression)
  models_glm_rdm7[ , paste0("class", i)] <- validation_group
  models_glm_rdm7[ , paste0("pred", i)] <- predict_validation
}

accrcy_glm_rdm7 <- as.matrix(accrcy_glm_rdm7)
accrcy_glm_rdm7m <- mean(accrcy_glm_rdm7)
models_glm_rdm7 <- models_glm_rdm7[,-1] # Remove the NaN column

#################################################################
write.csv(accrcy_glm_rdm7,paste(mainfolderpath, "glm_CV_validations_accry_bot7.csv", sep=""))
write.csv(models_glm_rdm7,paste(mainfolderpath, "glm_CV_validations_models_bot7.csv", sep=""))

#write.csv(accrcy_glm_rdm7,paste(mainfolderpath, "glm_CV_validations_accry_botbot7.csv", sep=""))      ## un-comment as necessary
#write.csv(models_glm_rdm7,paste(mainfolderpath, "glm_CV_validations_models_botbot7.csv", sep=""))     ## un-comment as necessary
#################################################################


#############################################
## Plot accuracy rates of different models ##
#############################################

## Load data
models_fi <- paste(mainfolderpath, "RandomForest_CV_validations_top50.csv", sep="")                             ## file path
models_dat <- read.table(models_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)            ## load to table
models_dat[,1] <- NULL

accry_fi <- paste(mainfolderpath, "RandomForest_CV_validations_accry_top50.csv", sep="")                        ## file path
accry_dat_top50 <- read.table(accry_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)        ## load to table
accry_dat_top50[,1] <- NULL

accry_fi <- paste(mainfolderpath, "RandomForest_CV_validations_accry_top7.csv", sep="")                         ## file path
accry_dat_top7 <- read.table(accry_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)         ## load to table
accry_dat_top7[,1] <- NULL

accry_dat_new <- accry_dat_top50
accry_dat_new[2,] <- accry_dat_top7
rownames(accry_dat_new) <- c("Top50","Top7") 
accry_dat_new <- data.frame(t(accry_dat_new))

accry_fi <- paste(mainfolderpath, "glm_CV_validations_accry_top7.csv", sep="")                                         ## file path
accry_dat_top7_glm <- read.table(accry_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)            ## load to table
accry_dat_top7_glm[,1] <- NULL
accry_dat_top7_glm <- data.frame(t(accry_dat_top7_glm))
colnames(accry_dat_top7_glm) <- c("Top7_glm")

accry_fi <- paste(mainfolderpath, "glm_CV_validations_accry_bot7.csv", sep="")                                         ## file path
accry_dat_bot7_glm <- read.table(accry_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)            ## load to table
accry_dat_bot7_glm[,1] <- NULL
accry_dat_bot7_glm <- data.frame(t(accry_dat_bot7_glm))
colnames(accry_dat_bot7_glm) <- c("Bot7_glm")

accry_fi <- paste(mainfolderpath, "glm_CV_validations_accry_botbot7.csv", sep="")                                         ## file path
accry_dat_botbot7_glm <- read.table(accry_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)            ## load to table
accry_dat_botbot7_glm[,1] <- NULL
accry_dat_botbot7_glm <- data.frame(t(accry_dat_botbot7_glm))
colnames(accry_dat_botbot7_glm) <- c("Botbot7_glm")

myData <- accry_dat_new %>% pivot_longer(everything())
head(myData)

ggplot(myData, aes(y= value, x = name, fill = name)) + 
  geom_boxplot(notch = FALSE, outlier.shape = NA,outlier.alpha = 0.1) + geom_jitter(width = 0.2) +
  theme_classic() + ylim(50, 100)

myData2 <- accry_dat_top7_glm %>% pivot_longer(everything())
myData3 <- rbind(myData, myData2)

ggplot(myData3, aes(y= value, x = name, fill = name)) + 
  geom_boxplot(notch = FALSE, outlier.shape = NA,outlier.alpha = 0.1) + geom_jitter(width = 0.2) +
  theme_classic() + ylim(25, 105)

myData4 <- data.frame("name"=NaN,"value"=NaN)
myData5 <- data.frame("value"=NaN)
for(j in 1:10) {
  accry_fi <- paste(mainfolderpath, "glm_CV_validations_accry_rdm7_", j, ".csv", sep="")
  accry_dat_rdm7_glm <- read.table(accry_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)         
  accry_dat_rdm7_glm[,1] <- NULL
  accry_dat_rdm7_glm <- data.frame(t(accry_dat_rdm7_glm))
  colnames(accry_dat_rdm7_glm) <- c("value")
  myData5 <- rbind(myData5, accry_dat_rdm7_glm)
  colnames(accry_dat_rdm7_glm) <- c(paste("Rdm7_glm", j, sep=""))
  dataT <- accry_dat_rdm7_glm %>% pivot_longer(everything())
  myData4 <- rbind(myData4, dataT)
}

myData4 <- myData4[-c(1), ]
myData5$name <- "Total"
myData5 <- (myData5[-c(1), ])
myData6 <- rbind(myData4, myData5)

ggplot(myData4, aes(y= value, x = name, fill = name)) + 
  geom_boxplot(notch = FALSE, outlier.shape = NA,outlier.alpha = 0.1) + geom_jitter(width = 0.2) +
  theme_classic() + ylim(50, 105)

ggplot(myData6, aes(y= value, x = name, fill = name)) + 
  geom_boxplot(notch = FALSE, outlier.shape = "o", outlier.alpha = 0.01) + geom_jitter(width = 0.1) +
  theme_classic() + ylim(50, 105) + theme(legend.position = "none")

myData_bot7 <- accry_dat_bot7_glm %>% pivot_longer(everything())
myData_botbot7 <- accry_dat_botbot7_glm %>% pivot_longer(everything())
myData7 <- rbind(myData4, myData_bot7,myData_botbot7)

ggplot(myData7, aes(y= value, x = name, fill = name)) + 
  geom_boxplot(notch = FALSE, outlier.shape = NA,outlier.alpha = 0.1) + geom_jitter(width = 0.2) +
  theme_classic() + ylim(25, 105)

myData8 <- rbind(myData7, myData3)

ggplot(myData8, aes(y= value, x = name, fill = name)) + 
  geom_boxplot(notch = FALSE, outlier.shape = NA,outlier.alpha = 0.1) + geom_jitter(width = 0.2) +
  theme_classic() + ylim(25, 105)

tiff(paste(mainfolderpath, "RF_glm_boxplots_1.tiff", sep=""), units="in", width=4, height=3, res=600)
ggplot(myData3, aes(y= value, x = name, fill = name)) + 
  geom_boxplot(notch = FALSE, outlier.shape = NA,outlier.alpha = 0.1) + geom_jitter(width = 0.2) +
  theme_classic() + ylim(25, 105)
dev.off()

tiff(paste(mainfolderpath, "RF_glm_boxplots_2.tiff", sep=""), units="in", width=16, height=3, res=600)
ggplot(myData7, aes(y= value, x = name, fill = name)) + 
  geom_boxplot(notch = FALSE, outlier.shape = NA,outlier.alpha = 0.1) + geom_jitter(width = 0.2) +
  theme_classic() + ylim(25, 105)
dev.off()


####################################################################
##  Correlation plots with clinical data included in tumor cases  ##
####################################################################

## Load data again
ageAdj_fi <- paste(mainfolderpath, "BiomarkerAnalysis_ageAdj.csv", sep="")                                      ## file path
ageAdj_dat <- read.table(ageAdj_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)            ## load to table
rownames(ageAdj_dat) <- ageAdj_dat[,1]                                                                          ## set rownames
ageAdj_dat[,1] <- NULL

## Separate clinical data and measurement data
clinical_dat <- (ageAdj_dat[,100:105])
ageAdj_dat <- ageAdj_dat[,-(92:105)]
ageAdj_dat <-t(ageAdj_dat)

## Lists taken from wilcox test data
wilcox_fi <- paste(mainfolderpath, "wilcox_test_results.csv", sep="")                                           ## file path
wilcox_dat <- read.table(wilcox_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)            ## load to table
wilcox_dat <- wilcox_dat[ , !(names(wilcox_dat) %in% c(""))]                                                    ## delete first column of numbers
wilcox_dat <- wilcox_dat[order(wilcox_dat$pAdj_w),]                                                             ## order based on p-values

## Top 20/50/70 protein markers
top20_proteins <- c(wilcox_dat$Pro[1:20])
top50_proteins <- c(wilcox_dat$Pro[1:50])
top70_proteins <- c(wilcox_dat$Pro[1:70])
top20_values <- ageAdj_dat[rownames(ageAdj_dat) %in% top20_proteins,]
top50_values <- ageAdj_dat[rownames(ageAdj_dat) %in% top50_proteins,]
top70_values <- ageAdj_dat[rownames(ageAdj_dat) %in% top70_proteins,]
top20_values <- t(top20_values)
top50_values <- t(top50_values)
top70_values <- t(top70_values)

## plot colours
col <- colorRampPalette(c("#4477AA","#77AADD","#FFFFFF","#EE9988","#BB4444"))

## Match ageAdjusted with clinical BY caseID
top20_mergedClinical <- merge(top20_values, clinical_dat, by = "row.names")
rownames(top20_mergedClinical) <- top20_mergedClinical[,1]; top20_mergedClinical[,1] <- NULL
top50_mergedClinical <- merge(top50_values, clinical_dat, by = "row.names")
rownames(top50_mergedClinical) <- top50_mergedClinical[,1]; top50_mergedClinical[,1] <- NULL
top70_mergedClinical <- merge(top70_values, clinical_dat, by = "row.names")
rownames(top70_mergedClinical) <- top70_mergedClinical[,1]; top70_mergedClinical[,1] <- NULL

top20_mergedClinical <- top20_mergedClinical[112:245,] # choose disease cases only (because Ctrl cases are NA for these variables)
top50_mergedClinical <- top50_mergedClinical[112:245,] # choose disease cases only (because Ctrl cases are NA for these variables)
top70_mergedClinical <- top70_mergedClinical[112:245,] # choose disease cases only (because Ctrl cases are NA for these variables)

## correlation matrices
top20_corrMatrix <- cor(top20_mergedClinical)
top50_corrMatrix <- cor(top50_mergedClinical)
top70_corrMatrix <- cor(top70_mergedClinical)

## Raw plots
top20_corrPlot <- corrplot(top20_corrMatrix, type="lower", col=col(200))
top50_corrPlot <- corrplot(top50_corrMatrix, type="lower", col=col(200))
top70_corrPlot <- corrplot(top70_corrMatrix, type="lower", col=col(200))

## Plots with hclust ordering
top20_corrPlotHclust <- corrplot(top20_corrMatrix, type="lower", order="hclust", col=col(200))
top50_corrPlotHclust <- corrplot(top50_corrMatrix, type="lower", order="hclust", col=col(200))
top70_corrPlotHclust <- corrplot(top70_corrMatrix, type="lower", order="hclust", col=col(200))

## Function for mapping p-values to matrix
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

## corr + signif matrices
top20_significance <- cor.mtest(top20_mergedClinical)
top50_significance <- cor.mtest(top50_mergedClinical)
top70_significance <- cor.mtest(top70_mergedClinical)

## Plot with all info
top20_signifHclust <- corrplot(top20_corrMatrix, col=col(200), order="hclust",
                               tl.col="black", tl.srt=45, diag=FALSE)
top50_signifHclust <- corrplot(top50_corrMatrix, col=col(200), order="hclust",
                               tl.col="black", tl.srt=45, diag=FALSE)
top70_signifHclust <- corrplot(top70_corrMatrix, col=col(200), order="hclust",
                               tl.col="black", tl.srt=45, diag=FALSE)

top20_corrMatrix2 <- top20_corrMatrix[1:20,1:20]
top50_corrMatrix2 <- top50_corrMatrix[1:50,1:50]
top70_corrMatrix2 <- top70_corrMatrix[1:70,1:70]

## Plot with all info
top20_signifHclust2 <- corrplot(top20_corrMatrix2, col=col(200), order="hclust",
                               tl.col="black", tl.srt=45, diag=FALSE)
top50_signifHclust2 <- corrplot(top50_corrMatrix2, col=col(200), order="hclust",
                               tl.col="black", tl.srt=45, diag=FALSE)
top70_signifHclust2 <- corrplot(top70_corrMatrix2, col=col(200), order="hclust",
                               tl.col="black", tl.srt=90, diag=FALSE, tl.cex = 0.6)

## save to CSV
write.csv(top20_corrMatrix, paste(mainfolderpath, "top20_corrMatrix.csv", sep=""))
write.csv(top50_corrMatrix, paste(mainfolderpath, "top50_corrMatrix.csv", sep=""))
write.csv(top70_corrMatrix, paste(mainfolderpath, "top70_corrMatrix.csv", sep=""))

## save to CSV
write.csv(top20_significance, paste(mainfolderpath, "top20_corrSignificance.csv", sep=""))
write.csv(top50_significance, paste(mainfolderpath, "top50_corrSignificance.csv", sep=""))
write.csv(top70_significance, paste(mainfolderpath, "top70_corrSignificance.csv", sep=""))

## spearman correlation matrices
top70_corrMatrix_S <- cor(top70_mergedClinical, method="spearman")
top70_corrMatrix2_S <- top70_corrMatrix_S[1:70,1:70]
top70_signifHclust2_S <- corrplot(top70_corrMatrix2_S, col=col(200), order="hclust",
                                tl.col="black", tl.srt=90, diag=FALSE, tl.cex = 0.6)

top50_corrMatrix_S <- cor(top50_mergedClinical, method="spearman")
top50_corrMatrix2_S <- top50_corrMatrix_S[1:50,1:50]
top50_signifHclust2_S <- corrplot(top50_corrMatrix2_S, col=col(200), order="hclust",
                                  tl.col="black", tl.srt=90, diag=FALSE, tl.cex = 0.6)

tiff(paste(mainfolderpath, "top50_corrMatrix_Sp.tiff", sep=""), units="in", width=12, height=12, res=600)
top50_signifHclust2_S <- corrplot(top50_corrMatrix2_S, col=col(200), order="hclust",
                                  tl.col="black", tl.srt=90, diag=FALSE, tl.cex = 1.0)
dev.off()


## Function for mapping p-values to matrix
cor.mtestSP <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], method="spearman", ...) ## kendall  spearman
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

## Top 80 protein markers
top80_proteins <- c(wilcox_dat$Pro[1:80])
top80_values <- ageAdj_dat[rownames(ageAdj_dat) %in% top80_proteins,]
top80_values <- t(top80_values)

## Match ageAdjusted with clinical by caseID
top80_mergedClinical <- merge(top80_values, clinical_dat, by = "row.names")
rownames(top80_mergedClinical) <- top80_mergedClinical[,1]; top80_mergedClinical[,1] <- NULL
top80_mergedClinical <- top80_mergedClinical[112:245,] # choose disease cases only (because Ctrl cases are NA for these variables)

top80_corrMatrix_Sp <- cor(top80_mergedClinical, method="spearman")
top80_significance_Sp_s <- cor.mtestSP(top80_mergedClinical)
pAdj_top80_significance_Sp_s <- matrix(p.adjust(as.vector(as.matrix(top80_significance_Sp_s)), method='fdr'),ncol=86)
rownames(pAdj_top80_significance_Sp_s) <- rownames(top80_corrMatrix_Sp)
colnames(pAdj_top80_significance_Sp_s) <- colnames(top80_corrMatrix_Sp)

write.csv(top80_corrMatrix_Sp, paste(mainfolderpath, "top80_corrMatrix_Sp.csv", sep=""))
write.csv(pAdj_top80_significance_Sp_s, paste(mainfolderpath, "top80_SPsignificance_pAdj.csv", sep=""))

n <- ncol(top80_corrMatrix_Sp)
SPsigList <- matrix(NA, (((n*n)-n)/2), 4)
cotr <- 1
for (i in 1:(80)) {
  for (j in (84):n) {
    if(pAdj_top80_significance_Sp_s[i,j]<0.05){
      SPsigList[cotr, 1] <- rownames(top80_corrMatrix_Sp)[i]
      SPsigList[cotr, 2] <- colnames(top80_corrMatrix_Sp)[j]
      SPsigList[cotr, 3] <- top80_corrMatrix_Sp[i,j]
      SPsigList[cotr, 4] <- pAdj_top80_significance_Sp_s[i,j]
      cotr <- cotr + 1
    }
  }
}
SPsigList <- SPsigList[-which(is.na(SPsigList[ ,1])), ]

write.csv(SPsigList, paste(mainfolderpath, "top80_SPsigList_short.csv", sep=""))


###########################################
##  Correlations between PEA and WB data ##
###########################################

## Load data
ageAdj_fi <- paste(mainfolderpath, "BiomarkerAnalysis_ageAdj.csv", sep="")                                      ## file path
ageAdj_dat <- read.table(ageAdj_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)            ## load to table
rownames(ageAdj_dat) <- ageAdj_dat[,1]                                                                          ## set rownames
ageAdj_dat[,1] <- NULL
ageAdj_values <- ageAdj_dat[1:134,1:(length(ageAdj_dat)-14)]

prot_fi <- paste(mainfolderpath, "ProteinAssayData.csv", sep="")
prot_dat <- read.csv(prot_fi, sep = ",", head = TRUE)
rownames(prot_dat) <- prot_dat$SampleID
prot_dat <- prot_dat[,-1]
prot_values <- prot_dat[1:(length(prot_dat)-7)]

# List of comparisons
tblcols <- expand.grid(1:ncol(prot_values), 1:ncol(ageAdj_values))

# Func for estimates
cfuncE <- function(var1, var2) {
  cor.test(prot_values[,var1], ageAdj_values[,var2], method="spearman",use="complete.obs")$estimate
}
# Func for p-values
cfuncP <- function(var1, var2) {
  cor.test(prot_values[,var1], ageAdj_values[,var2], method="spearman",use="complete.obs")$p.value
}
# Matrix for estimates
mat_est <- matrix(mapply(cfuncE, tblcols$Var1, tblcols$Var2), 
                  ncol = ncol(ageAdj_values), nrow = ncol(prot_values),
                  dimnames = list(colnames(prot_values), colnames(ageAdj_values)))
# Matrix for p-values
mat_pvals <- matrix(mapply(cfuncP, tblcols$Var1, tblcols$Var2), 
                    ncol = ncol(ageAdj_values), nrow = ncol(prot_values),
                    dimnames = list(colnames(prot_values), colnames(ageAdj_values)))

sig_HIF1a <- mat_pvals[,mat_pvals[1,] < 0.05]
sig_HIF2a <- mat_pvals[,mat_pvals[2,] < 0.05]
sig_pVHL <- mat_pvals[,mat_pvals[3,] < 0.05]
sig_Alk5FL <- mat_pvals[,mat_pvals[4,] < 0.05]
sig_Alk5ICD <- mat_pvals[,mat_pvals[5,] < 0.05]
sig_pSmad23 <- mat_pvals[,mat_pvals[6,] < 0.05]

# Correlations plot for WB readings vs TGFb-pVHL pathway proteins
list_fi <- paste(mainfolderpath, "WBvsPEAcorr_lists.csv", sep="")
list_dat <- read.csv(list_fi, sep = ",", head = TRUE)

mat_2plt <- data.frame(matrix(NaN, nrow(mat_est), nrow(list_dat["PEA_list"])))
colnames(mat_2plt) <- (list_dat[,2])
rownames(mat_2plt) <- rownames(mat_est)
for (i in 1:ncol(mat_2plt)){
  mat_2plt[1:6,i] <- (mat_est[,list_dat[,2][i]])
}

matp_2plt <- data.frame(matrix(NaN, nrow(mat_est), nrow(list_dat["PEA_list"])))
colnames(matp_2plt) <- (list_dat[,2])
rownames(matp_2plt) <- rownames(mat_est)

for (i in 1:ncol(matp_2plt)){
  matp_2plt[1:6,i] <- (mat_pvals[,list_dat[,2][i]])
}

col <- colorRampPalette(c("#2b6bd9","#4a78c7","#77AADD","#FFFFFF","#EE9988","#BB4444","#bb0c00"))

tiff(paste(mainfolderpath, "WB_vs_PEA_corrplot.tiff", sep=""), units="in", width=15, height=3, res=600)
wbvsPEAproteins <- corrplot(as.matrix(mat_2plt), method="circle", tl.col="black", tl.srt=90, tl.cex = 0.8,is.corr = FALSE,col.lim = c(-0.4,0.4),
                            p.mat = as.matrix(matp_2plt), insig = "label_sig", sig.level = c(.05), pch.cex = 2.0, pch.col = "gray20", pch='*',col=col(201))

dev.off()

matp_2plt_st <- ifelse(matp_2plt<0.05, "*", "")

tiff(paste(mainfolderpath, "WB_vs_PEA_heatmap.tiff", sep=""), units="in", width=15, height=3, res=600)
plot <- pheatmap(as.matrix(mat_2plt),
                 show_rownames = T, show_colnames = T, cluster_rows = T, cluster_cols = T, clustering_distance_cols = "euclidean",
                 trace = "none", density.info = "none", treeheight_row=0, treeheight_col=0, display_numbers = matp_2plt_st, fontsize_number=18,
                 cellheight=15, cellwidth=15, fontsize=15, width=20, height=20,angle_col=90, breaks=seq(-0.3,0.3,length.out=(100 + 1)))

plot
dev.off()

write.csv(mat_2plt, paste(mainfolderpath, "WB_vs_PEA_corr.csv", sep=""))
write.csv(matp_2plt, paste(mainfolderpath, "WB_vs_PEA_pvals.csv", sep=""))


## The end