%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create personalized liver model for Genotypes for contextualizing with gene
% expression data for each subject
% Partho Sen: 01.08.2020
% Updated: 01.08.2020

% 'HMR_0031': 0.19 1,2-diacylglycerol-LD-TAG pool[c] + 0.0014 1-acyl-PE pool[c] + 0.0024 1-radyl-2-acyl-sn-glycero-3-phosphocholine[c] + 0.0006 2-lysolecithin pool[c] + 0.0008 O-1-alk-1-enyl-2-acyl-sn-glycero-3-phosphoethanolamine[c] + 0.0092 PC-LD pool[c] + 0.0034 PE-LD pool[c] + 0.0016 PI pool[c] + 0.0002 PS-LD pool[c] + 0.0004 SM pool[c] + 0.44 TAG-LD pool[c] + 0.34 cholesterol-ester pool[c] + 0.005 cholesterol[c] + 0.005 fatty acid-LD-TG1 pool[c] => lipid droplet[c]
% Lipid droplet production set as Objective function for E-FFLUX
% Tutorial Cobra: https://github.com/opencobra/COBRA.tutorials/blob/master/dataIntegration/extractionTranscriptomic/tutorial_extractionTranscriptomic.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
javaaddpath('/Users/pasase/VirtualBox VMs/VM_share/SANELA_WABI/To_Partho_Sanela_working_copy/data_new_results/OPLS_B_Vs_C_F/MatlabExcelMac/Archive/jxl.jar');
javaaddpath('/Users/pasase/VirtualBox VMs/VM_share/SANELA_WABI/To_Partho_Sanela_working_copy/data_new_results/OPLS_B_Vs_C_F/MatlabExcelMac/Archive/MXL.jar');
import mymxl.*;
import jxl.*; 

%changeCobraSolver ('glpk', 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uu=0; 
disp("Getting the expression data as got from Olivier corrected for gender and batch...");

if uu
 %%% Remember here the unit is log-CPM from Limma and Voom (See Original
 %%% paper) + batch and sex adjusted!
 D = readtable('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq_Genotype/RNAseq_n216_normalised_batch_sex_Formatted.xlsx');
 Genes = table2array(D(:,1)); Genes(1)=[];
 Conditions = D.Properties.VariableNames; Conditions(1)=[];
 expression_all = table2array(D(:,2:size(D,2)));
 expression_all = str2double(expression_all); 
 % expression_all = 2.^str2double(expression_all); %% If converted into RAW values 
 %save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq_Genotype/Genes','Genes')
 %save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq_Genotype/Conditions','Conditions')
 %save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq_Genotype/expression_all','expression_all')
disp('done!')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load gene expression datasets i.e. log-CPM data for contextualization :)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Load and set gene expression datasets previously loaded!....')

load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq_Genotype/Genes') ; 
load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq_Genotype/Conditions');
load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq_Genotype/expression_all'); %%% log-CPM from Limma and Voom

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get the gene expression for the Liver genotypes or genetic variants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Read the header file: Same as the expression file :-')
[D,H] = xlsread('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq_Genotype/Header_RNAseq_n216_normalised_batch_sex.xlsx');
Subject_ID = H(1,:); Subject_ID(1) = [];
Grps = H(2,:);

uu=1;
if uu
disp('Select only exclusive cases and avoid mixed-up such as HetA vs HetB for a single subject'), %% 106 subjects
D(D==2) = 50; %% homo
D(D==1) = 10; %% hetero
D(D==0) = 0; %% WT
D(isnan(D))=0;
temp = sum(D); 

disp('Select only the subjectes either homo or heterozygous for one genotype only PNPLA3 or TM6 or HSD ...')
idx0 = find(temp==50); idx1 = find(temp==10); idx2 = find(temp==0); %% Add WTs for all the groups
idx = [idx0,idx1,idx2]; %% 149 exclusive gene variants only for PNPLA3, TMS6 or HSD
%idx = [idx0,idx1]; %% Combine homo and hetero but keep wild type separate

idx_WT_All = find(sum(D(:,idx))==0); % D(:,idx(idx_WT_All)); %% 43 exclusively wildtypes

%idx0_pnpla3 = find(D(1,idx)==0); %% Wild Type
idx1_pnpla3 = find(D(1,idx)==10); %% hetero %% D(:,idx(idx1_pnpla3)) %% 47 exclusively
idx2_pnpla3 = find(D(1,idx)==50); %% homo %% 22

%idx0_TMS6 = find(D(2,idx)==0); %% Wild Type
idx1_TMS6 = find(D(2,idx)==10); %% hetero %% D(:,idx(idx1_TMS6)) %% 13 exclusively
idx2_TMS6 = find(D(2,idx)==50); %% homo %% D(:,idx(idx2_TMS6)) %%  only 1 exclusively 

%idx0_HSD = find(D(3,idx)==0); %% Wild Type
idx1_HSD = find(D(3,idx)==10); %% hetero %% D(:,idx(idx1_HSD)) %% 20 exclusively
idx2_HSD = find(D(3,idx)==50); %% homo %% D(:,idx(idx2_HSD)) %% 3 exclusively

% idx0_pnpla3 = find(D(1,:)==0); %% Wild Type
% idx1_pnpla3 = find(D(1,:)==1); %% hetero
% idx2_pnpla3 = find(D(1,:)==2); %% homo
% 
% idx0_TMS6 = find(D(2,:)==0); %% Wild Type
% idx1_TMS6 = find(D(2,:)==1); %% hetero
% idx2_TMS6 = find(D(2,:)==2); %% homo
% 
% idx0_HSD = find(D(3,:)==0); %% Wild Type
% idx1_HSD = find(D(3,:)==1); %% hetero
% idx2_HSD = find(D(3,:)==2); %% homo


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Write genetic variants gene expression data to a File; %% N.B: Use R to
%%% Calculate the differential expression analysis and come back to Matlab

PNPLA3_WT_genes = expression_all(:,idx(idx_WT_All)); PNPLA3_WT_genes1 = [Conditions(idx(idx_WT_All));num2cell(PNPLA3_WT_genes)]; %% Nos. 89
PNPLA3_het_genes = expression_all(:,idx(idx1_pnpla3)); PNPLA3_het_genes1 = [Conditions(idx(idx1_pnpla3));num2cell(PNPLA3_het_genes)]; %% Nos. 89
PNPLA3_hom_genes = expression_all(:,idx(idx2_pnpla3)); PNPLA3_hom_genes1 = [Conditions(idx(idx2_pnpla3));num2cell(PNPLA3_hom_genes)]; %% Nos. 42
PNPLA3_hom_het = [PNPLA3_hom_genes1,PNPLA3_het_genes1]; size(PNPLA3_hom_het) %% 32523 X 69 %% Combine Het and Homo for power
PNPLA3_wt = PNPLA3_WT_genes1; size(PNPLA3_wt) %% 32523 X 43 %%

TM6_WT_genes  = expression_all(:,idx(idx_WT_All));  TM6_WT_genes1 = [Conditions(idx(idx_WT_All));num2cell(TM6_WT_genes)]; %% Nos. 19
TM6_het_genes = expression_all(:,idx(idx1_TMS6)); TM6_het_genes1 = [Conditions(idx(idx1_TMS6));num2cell(TM6_het_genes)]; %% Nos. 21
TM6_hom_genes = expression_all(:,idx(idx2_TMS6)); TM6_hom_genes1 = [Conditions(idx(idx2_TMS6));num2cell(TM6_hom_genes)]; %% Nos. 11
TM6_hom_het = [TM6_hom_genes1,TM6_het_genes1]; size(TM6_hom_het); %% 32523 X 43 %% 
TM6_wt = TM6_WT_genes1; size(TM6_wt) %% 32523 X 43 %%

HSD_WT_genes  = expression_all(:,idx_WT_All);  HSD_WT_genes1 = [Conditions(idx_WT_All);num2cell(HSD_WT_genes)]; %% Nos. 33
HSD_het_genes = expression_all(:,idx1_HSD); HSD_het_genes1 = [Conditions(idx1_HSD);num2cell(HSD_het_genes)]; %% Nos. 28
HSD_hom_genes = expression_all(:,idx2_HSD); HSD_hom_genes1 = [Conditions(idx2_HSD);num2cell(HSD_hom_genes)]; %% Nos. 15
HSD_hom_het = [HSD_hom_genes1,HSD_het_genes1]; size(HSD_hom_het) %% %% 32523 X 23 %%
HSD_wt = HSD_WT_genes1; size(HSD_wt) %% %% 32523 X 43 %%


uu=0;
if uu
 disp('Adding sheets to excel .... it might take a while crashes ....')  
 disp('Saturates Java heap space.... so need to write in a separate file ...')

 A = [['Gene_ID';Genes],PNPLA3_hom_het];
 A1 = [['Gene_ID';Genes],PNPLA3_wt];

 B = [['Gene_ID';Genes],TM6_hom_het];
 B1 = [['Gene_ID';Genes],TM6_wt];
 
 C = [['Gene_ID';Genes],HSD_hom_het];
 C1 = [['Gene_ID';Genes],HSD_wt];

writetable(table(A),'/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq_Genotype/Genotype_Expression_Combined_PNPLA3_hom_het_genes_ver2.csv')
writetable(table(A1),'/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq_Genotype/Genotype_Expression_Combined_PNPLA3_wt_genes_ver2.csv')

writetable(table(B),'/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq_Genotype/Genotype_Expression_Combined_TM6_hom_het_genes_ver2.csv')
writetable(table(B1),'/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq_Genotype/Genotype_Expression_Combined_TM6_wt_genes_ver2.csv')

writetable(table(C),'/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq_Genotype/Genotype_Expression_Combined_HSD_hom_het_genes_ver2.csv')
writetable(table(C1),'/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq_Genotype/Genotype_Expression_Combined_HSD_wt_genes_ver2.csv')

% xlwrite('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq_Genotype/Genotype_Expression_Combined_ver2.xlsx',A,'PNPLA3_Comb_genes')
% xlwrite('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq_Genotype/Genotype_Expression_Combined_ver2.xlsx',B,'TM6_Comb_genes')
% xlwrite('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq_Genotype/Genotype_Expression_Combined_ver2.xlsx',C,'HSD_Comb_genes')

disp('Keep a note of expression data in each of these files...')
disp('You may perform differential expression analysis based on different conditions on these datasets.... get the Rscript for diff analysis previous analysis of these projects')
disp('These datasets may also be used for the Reporter Metabolite analysis')

end

disp('End of RNAseq part... Moving to the Model reconstruction and QC and QA analysis')
disp('Ideally we should compare the genotypes hom & het with the WTs... but due to lack of statistical power.. we are selecting exclusively het or hom of PNPLA3/TM6/HSD subjects but WT for other genes...')
disp('Note 7 normal HCs subjects are WTs for all these three conditions')
disp('Now we can compare Gene variant cases with the controls, similarly as homo or het for a particular variant vs WT');
disp('Incuding WTs are also important because they show the NFALDs irrespective of het or homo');
disp('After running this block... format and combine the excel files and run the binding R script; alternatively the R script can be run later')
disp('Necessary step: manually combine the CSV files (heap size donot allow excel) to a common excel to  run the binding R script')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load reference models for contrxtualization in this case HMR2.0 cobra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Load reference liver GEMs both in Cobra format ....')

load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/HMR2_Cobra.mat'); %% Object created in the Test scripts. Format needed for contextualization.
Ref_HMR2_Cobra = HMR2_Cobra; 

disp('Setting Lipid droplet formation as objective function for E-flux');
Ref_HMR2_Cobra = setParam(Ref_HMR2_Cobra,'obj','HMR_0031',1); 

disp('Deciding on gene expression threshold and have one see ... Systemic evaluation paper in Cell Reports and Leif Varemo paper...');

%%% check if correct objective function is assigned
%%% printRxnFormula(Ref_HMR2_Cobra,'rxnAbbrList','HMR_0031','metNameFlag',1) 

% hist(expression_all,50); Plot a density plot in R
%%% Expression_threshold >= 0; %% log-scale anything greater than that
%%% exp(0) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% *** PNPLA3 section *** %%%
uu=1;
if uu
disp('Make variant models for PNPLA3...')
Cond = PNPLA3_hom_het(1,:); PNPLA3_hom_het(1,:)=[]; PNPLA3_hom_het = cell2mat(PNPLA3_hom_het); %Genes
disp('Create condition-specific Genotype liver GEMs PNPLA3 (Het,Homo) gene variants ..):')

[c,r] = size(PNPLA3_hom_het);

for ii=1:r
    
    disp(strcat('PNPLA3: ',num2str(ii),':',Cond{ii}))
    
    Exp_tar = PNPLA3_hom_het(:,ii);
    model = Personalized_NAFLD_iMAT_Eflux(Ref_HMR2_Cobra,Genes,Exp_tar,'/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Personalized_Models_Genotypes/',strcat('PNPLA3_',Cond{ii}));
    model
    
    disp('done!')

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% *** TM6 section *** %%%
uu=1;
if uu
disp('Make variant models for TMS6...')
Cond = TM6_hom_het(1,:); TM6_hom_het(1,:)=[]; TM6_hom_het = cell2mat(TM6_hom_het); %Genes

[c,r] = size(TM6_hom_het);

for ii=1:r
    
    disp(strcat('TM6: ',num2str(ii),':',Cond{ii}))
    
    Exp_tar = TM6_hom_het(:,ii);
    model = Personalized_NAFLD_iMAT_Eflux(Ref_HMR2_Cobra,Genes,Exp_tar,'/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Personalized_Models_Genotypes/',strcat('TM6_',Cond{ii}));
    model
    
    disp('done!')

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% *** HSD section *** %%%
uu=1;
if uu
disp('Make variant models for HSD...')
Cond = HSD_hom_het(1,:); HSD_hom_het(1,:)=[]; HSD_hom_het = cell2mat(HSD_hom_het); %Genes
disp('Create condition-specific Genotype liver GEMs HSD (Het,Homo) gene variants ..):')

[c,r] = size(HSD_hom_het);

for ii=1:r
    
    disp(strcat('HSD: ',num2str(ii),':',Cond{ii}))
    
    Exp_tar = HSD_hom_het(:,ii);
    model = Personalized_NAFLD_iMAT_Eflux(Ref_HMR2_Cobra,Genes,Exp_tar,'/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Personalized_Models_Genotypes/',strcat('HSD_',Cond{ii}));
    model
    
    disp('done!')

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% *** WT section: COMMON Wild Type: the pattern of the sequences mostly occuring among the population *** %%%
uu=1;
if uu
disp('Make WT models for PNPLA3...')
Cond = PNPLA3_wt(1,:); PNPLA3_wt(1,:)=[]; PNPLA3_wt = cell2mat(PNPLA3_wt); %Genes
disp('Create condition-specific Genotype liver GEMs PNPLA3 (WT) gene variants ..):')

[c,r] = size(PNPLA3_wt);

for ii=1:r
    
    disp(strcat('WT: ',num2str(ii),':',Cond{ii}))
    
    Exp_tar = PNPLA3_wt(:,ii);
    model = Personalized_NAFLD_iMAT_Eflux(Ref_HMR2_Cobra,Genes,Exp_tar,'/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Personalized_Models_Genotypes/',strcat('WT_',Cond{ii}));
    model
    
    disp('done!')

end
end

% Total time per run is from IBM_CPLEX MILP (solveMILPCobra) + functional_Model(consistency check 'fastcc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *** The END ***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%