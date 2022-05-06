#################################################################
# Factor Analysis
# Checking the site-effect on transcriptomics and metabolomics data
# Partho: April 19, 2021
# Edited: 02.05.2021
################################################################
#library(plotly)
library(gplots)
library(eulerr)
library(vegan)
library(ape)
library(stringr)
library(gplots)
library(ggplot2)
#library(ropls)
library(nlme) ## gives lme & nlme
library(lme4) ## gives lmer
library(gplots)
library(ggplot2)
library(ade4)
library(corrplot)
library(reshape)
library(Hmisc) ## Also using for imputed fumction
#library(ChemometricsWithR)
#library(ChemometricsWithRData)
library(plyr)
library(car)
library(maptools)
library(rgeos)
library(plotrix)
library(factoextra)
library(magrittr)
library(robustbase)
library(MASS)
library(scater)
library(SingleCellExperiment)
library(knitr)
require(gridExtra)
library(scRNAseq)
library(readxl)

options(stringsAsFactors = FALSE)

#require("maanova")
#require('mvnormtest')
#require(graphics)

rm(list=ls())
#setwd('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_MetaHit_PFAS/Scripts')
setwd('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Results/')


#################################### Factor analysis: Transcriptomics ################################################
load('../Data/RNAseq_Genotype/RNAseq_n216_normalised_batch_sex_Formatted.rda')

DH = read_excel('../Data/RNAseq_Genotype/Header_RNAseq_n216_normalised_batch_sex.xlsx',progress = readxl_progress(), .name_repair = "unique")
DH.CN <-  colnames(DH);
DH = data.frame(DH); colnames(DH) <- DH.CN; DH = DH[,-1]; DH = DH[1,];

DF.new = as.data.frame(read_excel('../Data/Clinical/RNAseq_centre_sex_formatted.xlsx',progress = readxl_progress(), .name_repair = "unique"));
idx = match(names(DH), DF.new[,1])

print(cbind(DF.new[idx,1], names(DH) ))

colnames(DFF) <- names(DH); #DFF = t(apply(DFF,2,as.numeric));
Annot = DF.new[idx,3:5];

DFF = t(scale(as.matrix(DFF), scale=F, center = F));
rownames(Annot) <- rownames(DFF)


pdf("../Results/Figure_new/Factors_Genes.pdf", width = 12, height =8, paper = 'special')

umi <- SingleCellExperiment(assays = list(counts =  2^t(DFF) ), colData = Annot )
#umi <- normalize(umi,return_log=1); # exprs(umi);
umi <- logNormCounts(umi); # exprs(umi);

drop_metabs <- apply(exprs(umi), 1, function(x) {var(x) == 0})
umi <- umi[!drop_metabs, ];

vars <- getVarianceExplained(umi)
vars1 = vars[,which(colnames(vars) %in% c("Center", 		
                                          "Sex", 
                                          "NAFLD_Group" ))];
print(head(vars1))


g=plotExplanatoryVariables(vars1, 
                           theme_size = 12, 
                           nvars_to_plot = length(colnames(vars1)), 
                           min_marginal_r2 = 0.01)
print(g)

dev.off()

#################################### Factor analysis: Metabolomics ################################################

##** Read old RNAseq file **##
DF.RNAseq  = read.table('../Data/RNAseq/cohort_rnaseq_pheno.csv',sep='\t', header = TRUE,check.names=FALSE) 
DF.RNAseq  = DF.RNAseq[1:nrow(DF.RNAseq),]; ##Check rows
DF.RNAseq  = DF.RNAseq[,c(1,3)]; ## Sample SPCID

##** Read Bar code matching file: as a part of clinical data **##
BC_Map  = read.table('../Data/Clinical/Barcodes_Matching_Clinical_Data_PS.csv',sep='\t', header = TRUE,check.names=FALSE) #  121 samples

idp = match(DF.RNAseq$SPICID, BC_Map[,2])
DF.RNAseq = cbind(DF.RNAseq,BC_Map[idp,]); ## ID's mapped into this!

Cers.DF = read.csv('../Data/Cers_Vits_Tuulia/Epos_vitamins_ceramides_11_2020_formatted.csv',header = T, sep=',',stringsAsFactors = F,check.names = F)
Cers.DF = t(Cers.DF); colnames(Cers.DF) = Cers.DF[1,]; Cers.DF = Cers.DF[-1,]; rowN = rownames(Cers.DF)
Cers.DF = apply(Cers.DF,2,as.numeric); rownames(Cers.DF) <- rowN;
Cers.DF = log2(Cers.DF);
#Cers.DF = data.frame(scale(Cers.DF,center = T, scale = T));
Cers.DF = apply(Cers.DF,2,function(x) tapply(x, rownames(Cers.DF), mean))


DF.RNAseq = DF.RNAseq[-which(is.na(DF.RNAseq$Barcode)), ];

idx1 = match(DF.RNAseq$Barcode, rownames(Cers.DF))

print(cbind( DF.RNAseq$Barcode[!is.na(idx1)], rownames(Cers.DF)[idx1[!is.na(idx1)]] ))

Annot = DF.RNAseq[!is.na(idx1),];
Cers.DF = Cers.DF[idx1[!is.na(idx1)],];

print(cbind( Annot$Barcode, rownames(Cers.DF) )); ### 42 subjects metabolomics overlap
rownames(Cers.DF) <-  Annot$Sample;


DF.new = as.data.frame(read_excel('../Data/Clinical/RNAseq_centre_sex_formatted.xlsx',progress = readxl_progress(), .name_repair = "unique"));
idx = match(rownames(Cers.DF), DF.new[,1])

print(cbind(DF.new[idx,1], rownames(Cers.DF) ))

Annot = DF.new[idx,3:5];
Cers.DF = scale(as.matrix(Cers.DF), scale=F, center = F);
rownames(Annot) <- rownames(Cers.DF);

Cers.DF = Cers.DF[-which(is.na(Annot$Center)),];
Annot = Annot[-which(is.na(Annot$Center)),];

print(cbind(rownames(Annot), rownames(Cers.DF) ))

pdf("../Results/Figure_new/Factors_Metabolites.pdf", width = 12, height =8, paper = 'special')

umi <- SingleCellExperiment(assays = list(counts =  2^t(Cers.DF) ), colData = Annot)
#umi <- normalize(umi,return_log=1); # exprs(umi);
umi <- logNormCounts(umi); # exprs(umi);

drop_metabs <- apply(exprs(umi), 1, function(x) {var(x) == 0})
umi <- umi[!drop_metabs, ];

vars <- getVarianceExplained(umi)
vars1 = vars[,which(colnames(vars) %in% c("Center", 		
                                          "Sex", 
                                          "NAFLD_Group" ))];
print(head(vars1))


g=plotExplanatoryVariables(vars1, 
                           theme_size = 12, 
                           nvars_to_plot = length(colnames(vars1)), 
                           min_marginal_r2 = 0.01)
print(g)

dev.off()



######################################## *********************************************************************** ###########################

