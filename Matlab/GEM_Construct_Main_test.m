%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here read curated or draft HTimmR 
% Get the gene sets for each reactions
% Create a vector of gene expression for each reactions
% Use Init to contexualize the each Th17 model
% This module is Steatosis specific contextualized GEMs
% Partho Sen 15th October 2017
% Updated and revised Sep 10, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

javaaddpath('/Users/pasase/VirtualBox VMs/VM_share/SANELA_WABI/To_Partho_Sanela_working_copy/data_new_results/OPLS_B_Vs_C_F/MatlabExcelMac/Archive/jxl.jar');
javaaddpath('/Users/pasase/VirtualBox VMs/VM_share/SANELA_WABI/To_Partho_Sanela_working_copy/data_new_results/OPLS_B_Vs_C_F/MatlabExcelMac/Archive/MXL.jar');
import mymxl.*;
import jxl.*; 

% Read HTimmR curated models %% Add the gene list

%HTimmR = importExcelModel('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_GEM_Tcell/Human_Models/HTimmR/HTmmiR_draft_Edited.xlsx');
%save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_GEM_Tcell/Human_Models/HTimmR/HTimmR_edited','HTimmR')

%load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_GEM_Tcell/Human_Models/HTimmR/HTiMMER_draft_model_copy.mat');
%HTimmR=HTiMMER_draft_model_copy;

load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/HTmmiR_Manual_Full_model.mat');
HTimmR=HTmmiR_Manual_Full_model;
HTimmR.description = 'Steatosis GEM';

%%%%%%%%%% Get genes for each reactions %%%%%%%%%%%%%%%%%

GenMat = full(HTimmR.rxnGeneMat);
GenesInRxns ={};
for uu =1:length(GenMat)
    
    uu
    GenesInRxns{uu} = HTimmR.genes(find(GenMat(uu,:))); %% List of genes per reactions
    
end



%%%%%%%%%% Contextualize with Steato gene expression dataset %%%%%%%%%%%%%%%%%
%%%%%%%%%% Create the reaction gene vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yy=0;
if yy
    
[D,H]=xlsread('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Steatosis_DEG_Dseq0.05.xlsx');
Genes = H(:,1); Genes(1)=[]; %% All these genes are differentially expressed between Seatosis and Control

%%%%%%%%%% Get the Mapping files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uu=0;
if uu
[G1,H1] = xlsread('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_GEM_Tcell/GeneID_Mapping/GeneID_Mapping');
disp('Convert Entrez to Ensemble... Skip if it is already in Ensemble format')
EnsembleID = Genes; 
EntrezID = G1(:,11); %% corresponding entrez ID

idp=[];

for yy = 1:length(Genes)
    
    yy

    idx = find(ismember(EntrezID,Genes(yy)));
    
    if(length(idx) > 0)
        
        idp = [idp; EnsembleID(idx(1))];
    end
    
    if(length(idx) == 0)
        
        idp = [idp;{'--'}];
    end
    
    
end

Genes = idp;
end



idx1=[];

for gg =1:length(GenesInRxns)
    
    gg
        idu = find(ismember(Genes,GenesInRxns{gg})); %% This could be improved for multiple genes
        
        if(length(idu) > 0)

            %idx1(gg) = mean(T_data(idu));
            idx1(gg) = 10;
      
        end
        
        if(length(idu) == 0)

            idx1(gg) = 0.0005;
            %idx1(gg) = 2;

      
        end
    
     
end

%idx1(isnan(idx1)) = -12; %% negative values
Steato_Rxns_Score = idx1;
  save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Steato_Rxns_Score','Steato_Rxns_Score');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run tInit for Steato draft model
%%%%%%%%%%%%%%%s%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Steato_Rxns_Score','Steato_Rxns_Score'); %% Load the reaction scores obtained from Gene expressions

Steato_Rxns_Score(Steato_Rxns_Score == 0) = -10; %% with zero means gene expression reaction removed
Steato_Rxns_Score(isnan(Steato_Rxns_Score)) = -10; %% negative values

%%%% tINIT: Method I %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% getInit (upper level functions) and runInit may be together called tINIT

[D,H]=xlsread('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Results/MetabolitesScores_4_GEMs_formatted.xlsx','Measured_Metabolites');
Met_Comp = H(:,1); Met_Comp(1)=[];
Met_Keep_Use = H(:,2); Met_Keep_Use(1)=[];
Met_Keep_Use = unique(Met_Keep_Use); %%% Presence of metabolites added as a constrain to the model 

HTimmR.rxnNames = HTimmR.rxns;

tt=1; %% Amazing it works! Also with Cobra model struct..
if tt
    tic
          setRavenSolver('mosek'); %% look RAVEN documentation to execute this !
          
          [HTimmR_Simp,deletedReactions, deletedMetabolites]=simplifyModel(HTimmR); %% remove exchange reaction for contexualization
          [Steato_draft_model,Steato_deletedRxns,Steato_metProduction,Steato_fValue] = runINIT(HTimmR_Simp,Steato_Rxns_Score',Met_Keep_Use);
          
          Sol=solveLP(Steato_draft_model,'max');

          save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Steato_draft_model','Steato_draft_model');
          save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Steato_deletedRxns','Steato_deletedRxns');

          exportToExcelFormat(Steato_draft_model,'/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Steato_draft_model.xlsx')

     toc
     
end



return
