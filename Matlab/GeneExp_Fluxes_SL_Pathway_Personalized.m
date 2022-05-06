%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gene and Fluxes GSL pathway
% Partho Sen: 08.08.2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

javaaddpath('/Users/pasase/VirtualBox VMs/VM_share/SANELA_WABI/To_Partho_Sanela_working_copy/data_new_results/OPLS_B_Vs_C_F/MatlabExcelMac/Archive/jxl.jar');
javaaddpath('/Users/pasase/VirtualBox VMs/VM_share/SANELA_WABI/To_Partho_Sanela_working_copy/data_new_results/OPLS_B_Vs_C_F/MatlabExcelMac/Archive/MXL.jar');
import mymxl.*;
import jxl.*; 

[E1,G1] = xlsread('../Data/RNAseq_NAFLD/NAFLD_differential_expression_DSeq2_ver2.0.xlsx','NASH(F0-1) vs. NAFL'); G1 = G1(2:length(G1(:,1)),1);
[E2,G2] = xlsread('../Data/RNAseq_NAFLD/NAFLD_differential_expression_DSeq2_ver2.0.xlsx','NASH(F2) vs. NAFL'); G2 = G2(2:length(G2(:,1)),1);
[E3,G3] = xlsread('../Data/RNAseq_NAFLD/NAFLD_differential_expression_DSeq2_ver2.0.xlsx','NASH(F3) vs. NAFL'); G3 = G3(2:length(G3(:,1)),1);
[E4,G4] = xlsread('../Data/RNAseq_NAFLD/NAFLD_differential_expression_DSeq2_ver2.0.xlsx','NASH(F4) vs. NAFL'); G4 = G4(2:length(G4(:,1)),1);

[E5,G5] = xlsread('../Data/RNAseq_NAFLD/NAFLD_differential_expression_DSeq2_ver2.0.xlsx','NASH(F2-3) vs. NASH(F0-1)'); G5 = G5(2:length(G5(:,1)),1);
[E6,G6] = xlsread('../Data/RNAseq_NAFLD/NAFLD_differential_expression_DSeq2_ver2.0.xlsx','NASH(F3-4) vs. NAFL+NASH(F0-1)'); G6 = G6(2:length(G6(:,1)),1);
[E7,G7] = xlsread('../Data/RNAseq_NAFLD/NAFLD_differential_expression_DSeq2_ver2.0.xlsx','NASH(F3-4) vs. NAFL+NASH(F0-2)'); G7 = G7(2:length(G7(:,1)),1);

G  = unique([G1;G2;G3;G4;G5;G6;G7]);
Mat = zeros(length(G),7);
Mat_Pval = zeros(length(G),7);

[tf,loc] = ismember(G,G1); [~,p] = sort(loc(tf)); idx = find(tf); idx1 = idx(p);
[tf,loc] = ismember(G,G2); [~,p] = sort(loc(tf)); idx = find(tf); idx2 = idx(p);
[tf,loc] = ismember(G,G3); [~,p] = sort(loc(tf)); idx = find(tf); idx3 = idx(p);
[tf,loc] = ismember(G,G4); [~,p] = sort(loc(tf)); idx = find(tf); idx4 = idx(p);
[tf,loc] = ismember(G,G5); [~,p] = sort(loc(tf)); idx = find(tf); idx5 = idx(p);
[tf,loc] = ismember(G,G6); [~,p] = sort(loc(tf)); idx = find(tf); idx6 = idx(p);
[tf,loc] = ismember(G,G7); [~,p] = sort(loc(tf)); idx = find(tf); idx7 = idx(p);

Mat(idx1,1) = E1(:,2);
Mat(idx2,2) = E2(:,2);
Mat(idx3,3) = E3(:,2);
Mat(idx4,4) = E4(:,2);
Mat(idx5,5) = E5(:,2);
Mat(idx6,6) = E6(:,2);
Mat(idx7,7) = E7(:,2);

Mat_Pval(idx1,1) = E1(:,6);
Mat_Pval(idx2,2) = E2(:,6);
Mat_Pval(idx3,3) = E3(:,6);
Mat_Pval(idx4,4) = E4(:,6);
Mat_Pval(idx5,5) = E5(:,6);
Mat_Pval(idx6,6) = E6(:,6);
Mat_Pval(idx7,7) = E7(:,6);

Mat = 2.^(Mat); %% logs to abs %%

load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/HMR2_RAVEN.mat'); %% Object created in the Test scripts. Format needed for contextualization.
model = HMR2_RAVEN; 

Rxn_Gidx = zeros(length(model.rxns),7);
Rxn_Gid_Genes = {};

for ii=1:length(model.rxns)
    ii
    
    g = model.genes(find(model.rxnGeneMat(ii,:)));
    [tf,loc] = ismember(G,g); [~,p] = sort(loc(tf)); idx = find(tf); idx0 = idx(p);

    if(length(idx0)>1)
          Rxn_Gidx(ii,:) = mean(Mat(idx0,:),1);
          Rxn_Gid_Genes{ii} = G(idx0);
    end
    if(length(idx0)==1)
        
          Rxn_Gidx(ii,:) = Mat(idx0,:);
          Rxn_Gid_Genes{ii} = G(idx0);

    end

end

Subs = ["Sphingolipid metabolism","Glycosphingolipid metabolism","Glycosphingolipid biosynthesis-globo series","Glycosphingolipid biosynthesis-ganglio series","Glycosphingolipid biosynthesis-lacto and neolacto series"];
[tf,loc] = ismember(string(model.subSystems),Subs); [~,p] = sort(loc(tf)); idu = find(tf); idu = idu(p);

labss = {'NASH(F01) vs. NAFL','NASH(F2) vs. NAFL','NASH(F3) vs. NAFL','NASH(F4) vs. NAFL',...
    'NASH(F2-3) vs. NAFL-NASH(F0-1)','NASH(F3-4) vs. NAFL-NASH(F0-1)','NASH(F3-4) vs. NAFL-NASH(F0-2)'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Gene expression %%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rxn_Genes = log2(Rxn_Gidx(idu,:)); Rxn_Genes(isinf(Rxn_Genes)) = 0;
Rnoms = model.rxns(idu);
Gns = Rxn_Gid_Genes(idu);

% filter %%
FC1 = Rxn_Genes(any(abs(Rxn_Genes)' > 0),:);
Rnoms2 = Rnoms(any(abs(Rxn_Genes)' > 0));
Eqs1 = constructEquations(model,Rnoms2,0);
Eqs0 = Eqs1;

[Eqs1,id]= unique(Eqs1,'stable');
FC1 = FC1(id,:);
Rnoms2 = Rnoms2(id);

figure(8789)
imagesc(FC1)

set(gca,'fontsize',6,'Ytick',[1:length(Eqs1)],'YTickLabel',Eqs1,'Xtick',[1:7],...
    'XTickLabel',labss,'fontweight','bold');
colormap(bluewhitered)
colorbar
xtickangle(45)

disp('Saving GL gene expression...')
FC11 = [['Gl_Rxns',labss(3:7)];[Eqs1,num2cell(FC1(:,3:7))]];
    %xlwrite('../Results/Simulations/GeneExp_Fluxes_GSLs_Personalized_NAFLD.xlsx',FC11,'GeneExp_FCs');
disp('done!')


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read flux states %%%%%%%%%%%%%%%%%%%%%%%%%%%%

uu=1; %% switch off and on
if uu
    
disp('Loading personalized NAFLD models into a structure....Please wait!')

resPath2 ='/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Personalized_Models/';
files = dir(resPath2);

disp('Getting personalized models from the directory...')

for K = 1:length(files)
  thisdir = files(K).name;
  subdirinfo{K} = thisdir;
end
subdirinfo(1:3) =[]; %% remove the unwanted

disp('Loading personalized models from the directory...')

fooNames = subdirinfo; 
modelsP_NAFLD={}; 

for pp=1:length(fooNames) 
    
    strcat('Model: ',num2str(pp))
    
u=strcat(resPath2,fooNames(pp));
load(char(u));
modelsP_NAFLD{pp,1} = model;

end

model_names = fooNames;
REF = modelsP_NAFLD;
end

uRxns =  unique(Rnoms2);
Flx = zeros(length(uRxns),length(REF));
G_flux = [];

disp('Performing simulations ....')

for jj=1:length(REF)
    jj
    
model = REF{jj};

[rxns_uptake, uptIDX] = getExchangeRxns(model,'in'); %% -1000 - 1000
[rxns_prod, prodIDX] = getExchangeRxns(model,'out'); %% 0 - 1000
[rxns_both, bothIDX] = getExchangeRxns(model,'both'); %% -1000 - 1000 %% All exchange reactions

model = setParam(model,'ub',prodIDX,10); %% Cannot uptake but can secrete

model = setParam(model,'ub', {'HMR_9034','HMR_9047','HMR_9048','HMR_9058','HMR_0031',...
    },[1000,1000,1000,1000,1000]); %  glucose uptake, water uptake/produce, oxygen uptake, CO2 produced (direction of CO2 is different)
model = setParam(model,'lb', {'HMR_9034','HMR_9047','HMR_9048','HMR_9058','HMR_0031',...
    },[0,-1000,0,0,0]); % 

model = setParam(model,'obj',{'HMR_0031'},1); %% Optimization of the biomass %%

[tf,loc] = ismember(uRxns,model.rxns);
%[model.rxns(loc(find(loc))),uRxns(find(loc))]

%uRxns1 = uRxns(find(loc));
[minFlux,Sol] = fastFVA(model,90,'max','ibm_cplex',uRxns(find(loc)));

Flx(find(loc),jj) = Sol;

end

%F = log2(Flx(find(sum(Flx,2)),:)); F(~isfinite(F)) = 0;
%F = Flx(find(sum(Flx,2)),:); 

Eqs2 = constructEquations(HMR2_RAVEN,uRxns);

disp('Saving Fluxes...')

FLX = [['Sl_Rxns',model_names];[Eqs2,num2cell(Flx)]];
  %xlwrite('../Results/Simulations/GeneExp_Fluxes_GSLs_Personalized_NAFLD.xlsx',FLX,'Fluxes_model');
disp('done!')
