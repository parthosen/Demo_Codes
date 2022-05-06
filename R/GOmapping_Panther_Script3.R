##########################################################################################
# EPoS project
# Partho Sen: Jan 13 2019
# RNAseq: This Script reads MetaData files and gets the Counts file to make a common Table
# 64253 transcripts standarized #
# T1=cbind(T1,0); T2=t(apply(T1, 1, function(x) tapply(x, colnames(T1), mean)))
# DEG both by DSEQ and ttest #
###########################################################################################
library('piano')
library('biomaRt')
library('gplots')
library('calibrate')
library('RColorBrewer')

rm(list=ls())
setwd('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Results/')

DFF  = read.table('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq/cohort_rnaseq_pheno.csv',sep='\t', header = TRUE,check.names=FALSE) 
DFF  = DFF [1:nrow(DFF),]; ##Check rows

Samp = DFF[,1];
FN = DFF[,2];
SPICD = DFF[,3];
Grp = DFF[,5];


NC2 = read.table("NASH_DEG_Dseq0.05.csv",sep=',',header = TRUE,check.names=FALSE)
SC2  = read.table("Steatosis_DEG_Dseq0.05.csv",sep=',',header = TRUE,check.names=FALSE)
CC2  = read.table("Cirrhosis_DEG_Dseq0.05.csv",sep=',', header = TRUE,check.names=FALSE)


################### ****** RUN Panther Classifications ************ ####################################

source('../Scripts/pather.plots.R')

### *** NASH Vs. Ctrl *** #####
pather.plots("./pantherChart_BP_NASH.txt","./Figure/pantherChart_BP_NASH_bar.pdf","./Figure/pantherChart_BP_NASH_pie.pdf","NASH Vs Ctrl","Biological process (BP)");
pather.plots("./pantherChart_PC_NASH.txt","./Figure/pantherChart_PC_NASH_bar.pdf","./Figure/pantherChart_PC_NASH_pie.pdf","NASH Vs Ctrl","Protein Classes");
pather.plots("./pantherChart_Pathways_NASH.txt","./Figure/pantherChart_Pathways_NASH_bar.pdf","./Figure/pantherChart_Pathways_NASH_pie.pdf","NASH Vs Ctrl","Human Pathways (non-specific)")

### *** Steatosis Vs. Ctrl *** #####
pather.plots("./pantherChart_BP_Steatosis.txt","./Figure/pantherChart_BP_Steatosis_bar.pdf","./Figure/pantherChart_BP_Steatosis_pie.pdf","Steatosis Vs Ctrl","Biological process (BP)");
pather.plots("./pantherChart_PC_Steatosis.txt","./Figure/pantherChart_PC_Steatosis_bar.pdf","./Figure/pantherChart_PC_Steatosis_pie.pdf","Steatosis Vs Ctrl","Protein Classes");
pather.plots("./pantherChart_Pathways_Steatosis.txt","./Figure/pantherChart_Pathways_Steatosis_bar.pdf","./Figure/pantherChart_Pathways_Steatosis_pie.pdf","Steatosis Vs Ctrl","Human Pathways (non-specific)")

### *** Cirrhosis Vs. Ctrl *** #####

pather.plots("./pantherChart_BP_Cirrhosis.txt","./Figure/pantherChart_BP_Cirrhosis_bar.pdf","./Figure/pantherChart_BP_Cirrhosis_pie.pdf","Cirrhosis Vs Ctrl","Biological process (BP)");
pather.plots("./pantherChart_PC_Cirrhosis.txt","./Figure/pantherChart_PC_Cirrhosis_bar.pdf","./Figure/pantherChart_PC_Cirrhosis_pie.pdf","Cirrhosis Vs Ctrl","Protein Classes");
pather.plots("./pantherChart_Pathways_Cirrhosis.txt","./Figure/pantherChart_Pathways_Cirrhosis_bar.pdf","./Figure/pantherChart_Pathways_Cirrhosis_pie.pdf","Cirrhosis Vs Ctrl","Human Pathways (non-specific)")







#################################################################### *************************** The END *********************** #############################################################################