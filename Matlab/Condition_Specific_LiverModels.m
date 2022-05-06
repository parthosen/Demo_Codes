%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contextualize liver models based on different conditions
% Partho Sen: 17th Jan 2020
% Updated and revised 17th Jan 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
javaaddpath('/Users/pasase/VirtualBox VMs/VM_share/SANELA_WABI/To_Partho_Sanela_working_copy/data_new_results/OPLS_B_Vs_C_F/MatlabExcelMac/Archive/jxl.jar');
javaaddpath('/Users/pasase/VirtualBox VMs/VM_share/SANELA_WABI/To_Partho_Sanela_working_copy/data_new_results/OPLS_B_Vs_C_F/MatlabExcelMac/Archive/MXL.jar');
import mymxl.*;
import jxl.*; 

% uu=0;
% if uu
% D = readtable('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq_Genotype/RNAseq_n216_normalised_batch_sex_Formatted.xlsx');
% Genes = table2array(D(:,1)); Genes(1)=[];
% Conditions = D.Properties.VariableNames; Conditions(1)=[];
% expression_all = table2array(D(:,2:size(D,2)));
% %expression_all = 2.^str2double(expression_all);
% expression_all = abs(str2double(expression_all)); %% absolute magnitude of the expression values
% 
% save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq_Genotype/Genes','Genes')
% save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq_Genotype/Conditions','Conditions')
% save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq_Genotype/expression_all','expression_all')
% end

disp('Load and set gene expression datasets ....')

load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq_Genotype/Genes')
load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq_Genotype/Conditions')
load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq_Genotype/expression_all')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Get the group/condition-wise indices for the gene expression/RNAseq....')

ind_normal    = regexp(Conditions, 'Normal');    ind_normal = find(~cellfun(@isempty,ind_normal));
ind_Steatosis = regexp(Conditions, 'Steatosis'); ind_Steatosis = find(~cellfun(@isempty,ind_Steatosis));
ind_NASH_F1 = regexp(Conditions, 'NASH_F0_1');   ind_NASH_F1 = find(~cellfun(@isempty,ind_NASH_F1));
ind_NASH_F2 = regexp(Conditions, 'NASH_F2');     ind_NASH_F2 = find(~cellfun(@isempty,ind_NASH_F2));
ind_NASH_F3 = regexp(Conditions, 'NASH_F3');     ind_NASH_F3 = find(~cellfun(@isempty,ind_NASH_F3));
ind_NASH_F4 = regexp(Conditions, 'Cirrhosis_F4'); ind_NASH_F4 = find(~cellfun(@isempty,ind_NASH_F4));
disp('done!')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load reference liver GEMs both in RAVEN and Cobra formats ....')

load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/HMR2_Cobra.mat'); %% Object created in the Test scripts. Format needed for contextualization.
load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/HMR2_RAVEN.mat'); %% Object created in the Test scripts. Format needed to Fit and Check tasks.

Ref_HMR2_Cobra = HMR2_Cobra; %% For contextualization using GIMME and E-Flux
Ref_HMR2_RAVEN = HMR2_RAVEN; %% For gapFilling, Check and Fit Taks; Model QCs and Curations

Ref_HMR2_RAVEN.id = "HMR_2.00_RAVEN";
Ref_HMR2_Cobra.id = "HMR_2.00_Cobra";

Ref_HMR2_RAVEN.c(find(ismember(Ref_HMR2_RAVEN.rxns,'biomass_components')))=1; %% may not be used here 
Ref_HMR2_Cobra.c(find(ismember(Ref_HMR2_Cobra.rxns,'biomass_components')))=1; %% set biomass components pre-requisite for the GIMME 
disp('done!')

%%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%   Contextualize and create a condition-sepcific model using the expression datasets %
%%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
uu=0;
if uu
disp('Contextualize HMR 2.00 to 6 different coniditons [NAFLD?Spectrum]');

%%% *** Healthy Controls *** %%%
disp('Create condition-specific liver GEMs for healthy controls (HC)/normals')
[model_hc,RxnScores_hc,parsedGPR_hc] = ConditionSpecificEfluxGIMME(Ref_HMR2_Cobra,Genes,expression_all,ind_normal,'HC_GIMME_Eflux');
save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/Cobra_no_QC/model_hc.mat','model_hc'); %% Cobra format 
disp('done!');

%%% *** Steatosis *** %%%
disp('Create liver GEMs for Steatosis')
[model_stea,RxnScores_stea,parsedGPR_stea] = ConditionSpecificEfluxGIMME(Ref_HMR2_Cobra,Genes,expression_all,ind_Steatosis,'Steatosis_GIMME_Eflux');
save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/Cobra_no_QC/model_stea.mat','model_stea'); %% Cobra format 
disp('done!');

%%% *** NASH_F0/1 *** %%%
disp('Create liver GEMs for NASH_F0/1')
[model_nashF01,RxnScores_nashF01,parsedGPR_nashF01] = ConditionSpecificEfluxGIMME(Ref_HMR2_Cobra,Genes,expression_all,ind_NASH_F1,'NASH_F01_GIMME_Eflux');
save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/Cobra_no_QC/model_nashF01.mat','model_nashF01'); %% Cobra format 
disp('done!');

%%% *** NASH_F2 *** %%%
disp('Create liver GEMs for NASH_F2')
[model_nashF2,RxnScores_nashF2,parsedGPR_nashF2] = ConditionSpecificEfluxGIMME(Ref_HMR2_Cobra,Genes,expression_all,ind_NASH_F2,'NASH_F2_GIMME_Eflux');
save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/Cobra_no_QC/model_nashF2.mat','model_nashF2'); %% Cobra format 
disp('done!');

%%% *** NASH_F3 *** %%%
disp('Create liver GEMs for NASH_F3')
[model_nashF3,RxnScores_nashF3,parsedGPR_nashF3] = ConditionSpecificEfluxGIMME(Ref_HMR2_Cobra,Genes,expression_all,ind_NASH_F3,'NASH_F3_GIMME_Eflux');
save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/Cobra_no_QC/model_nashF3.mat','model_nashF3'); %% Cobra format 
disp('done!');

%%% *** NASH_F4 (Cirrhosis) *** %%%
disp('Create liver GEMs for NASH_F4(Cirrhosis)')
[model_nashF4,RxnScores_nashF4,parsedGPR_nashF4] = ConditionSpecificEfluxGIMME(Ref_HMR2_Cobra,Genes,expression_all,ind_NASH_F4,'NASH_F4_GIMME_Eflux');
save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/Cobra_no_QC/model_nashF4.mat','model_nashF4'); %% Cobra format 
disp('done!');
end

%%%%%%%%%%%%%%%%%%%%% Load condition-specific cobra models and convert to RAVEN formats, for QC analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load Models (Cobra formats) and convert to RAVEN formats for the QC analysis')
uu=0;
if uu
load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/Cobra_no_QC/model_hc.mat'); model_hc = ravenCobraWrapper(model_hc); model_hc.id = 'HC_GIMME_Eflux';
load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/Cobra_no_QC/model_stea.mat'); model_stea = ravenCobraWrapper(model_stea); model_stea.id = 'Steatosis_GIMME_Eflux';
load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/Cobra_no_QC/model_nashF01.mat'); model_nashF01 = ravenCobraWrapper(model_nashF01); model_nashF01.id = 'NASH_F01_GIMME_Eflux';
load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/Cobra_no_QC/model_nashF2.mat'); model_nashF2 = ravenCobraWrapper(model_nashF2); model_nashF2.id = 'NASH_F2_GIMME_Eflux';
load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/Cobra_no_QC/model_nashF3.mat'); model_nashF3 = ravenCobraWrapper(model_nashF3); model_nashF3.id = 'NASH_F3_GIMME_Eflux';
load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/Cobra_no_QC/model_nashF4.mat'); model_nashF4 = ravenCobraWrapper(model_nashF4); model_nashF4.id = 'NASH_F4_GIMME_Eflux';

save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/model_hc.mat','model_hc');
save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/model_stea.mat','model_stea');
save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/model_nashF01.mat','model_nashF01');
save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/model_nashF2.mat','model_nashF2');
save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/model_nashF3.mat','model_nashF3');
save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/model_nashF4.mat','model_nashF4');
end

%%%%%%%%%%%%%%%%%%%%% Load condition-specific cobra models and convert to RAVEN formats, for QC analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%
uu=0;
if uu
disp('Load Models (now RAVEN formats) for the QC analysis....')

load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/model_hc.mat');
load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/model_stea.mat');
load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/model_nashF01.mat');
load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/model_nashF2.mat');
load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/model_nashF3.mat');
load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/model_nashF4.mat');

disp('Group the models into a cell array....')

REF{1} = model_hc; REF{2} = model_stea; REF{3} = model_nashF01; REF{4} = model_nashF2; REF{5} = model_nashF3; REF{6} = model_nashF4;

%%%%%%%%%%%%%%%%%%%%% Necessary Step: Gap filling and fit tasks section by comparing HMR2.00/HTimmR %%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Adding Essentials.... ')
ref_model{1} = HMR2_RAVEN; %% group reference model 'HMR_RAVEN' into a Cell array!
for ii=1:6
ii    
REF{ii}.metComps = ref_model{1}.metComps(ismember(ref_model{1}.mets,strtok(cellstr(REF{ii}.mets),'[')));
REF{ii}.metNames = ref_model{1}.metNames(ismember(ref_model{1}.mets,strtok(cellstr(REF{ii}.mets),'[')));
end


    
setRavenSolver('mosek');

disp("Fit tasks for healthy controls (HC) .... ");
[model_hc,addedRxns]=fitTasks(REF{1},ref_model{1},'/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Scripts/Metabolic_Task2.xlsx',true);
save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/model_hc.mat','model_hc');
exportToExcelFormat(model_hc,'/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/xls/model_hc.xlsx')

disp('done!')

disp("Fit tasks for steatosis .... ");
[model_stea,addedRxns]=fitTasks(REF{2},ref_model{1},'/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Scripts/Metabolic_Task2.xlsx',true);
save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/model_stea.mat','model_stea');
exportToExcelFormat(model_hc,'/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/xls/model_stea.xlsx')

disp('done!')

disp("Fit tasks for NASH_F0/1 .... ");
[model_nashF01,addedRxns]=fitTasks(REF{3},ref_model{1},'/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Scripts/Metabolic_Task2.xlsx',true);
save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/model_nashF01.mat','model_nashF01');
exportToExcelFormat(model_hc,'/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/xls/model_nashF01.xlsx')

disp('done!')

disp("Fit tasks for NASH_F2 .... ");
[model_nashF2,addedRxns]=fitTasks(REF{4},ref_model{1},'/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Scripts/Metabolic_Task2.xlsx',true);
save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/model_nashF2.mat','model_nashF2');
exportToExcelFormat(model_hc,'/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/xls/model_nashF2.xlsx')

disp('done!')

disp("Fit tasks for NASH_F3 .... ");
[model_nashF3,addedRxns]=fitTasks(REF{5},ref_model{1},'/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Scripts/Metabolic_Task2.xlsx',true);
save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/model_nashF3.mat','model_nashF3');
exportToExcelFormat(model_hc,'/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/xls/model_nashF3.xlsx')
disp('done!')

disp("Fit tasks for NASH_F4 .... ");
[model_nashF4,addedRxns]=fitTasks(REF{6},ref_model{1},'/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Scripts/Metabolic_Task2.xlsx',true);
save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/model_nashF4.mat','model_nashF4');
exportToExcelFormat(model_hc,'/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/xls/model_nashF4.xlsx')
disp('done!')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Manual Check: Simulation for each model %%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Load Models (now RAVEN formats) for Simulations ....')

load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/model_hc.mat');
load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/model_stea.mat');
load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/model_nashF01.mat');
load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/model_nashF2.mat');
load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/model_nashF3.mat');
load('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/RAVEN_QC/model_nashF4.mat');

REF{1} = model_hc; REF{2} = model_stea; REF{3} = model_nashF01; REF{4} = model_nashF2; REF{5} = model_nashF3; REF{6} = model_nashF4;

setRavenSolver('mosek')

prompt = 'Select a model for simulation!Enter a number between 1 to 6:';
x = input(prompt); %% Enter the model number or the conditions from here:-

model_ind = x; %% Chnage here for the conditions
model = REF{model_ind};  
% model = setParam(model,'ub',getExchangeRxns(model),0); %% Note model should not uptake or produce anything unless allowed!
% model = setParam(model,'lb',getExchangeRxns(model),0);

model = setParam(model,'ub', {'HMR_9034','HMR_9047','HMR_9048','HMR_9058'},[100,1000,1000,1000]); % glucose, water, oxygen
model = setParam(model,'lb', {'HMR_9034','HMR_9047','HMR_9048','HMR_9058'},[-100,-1000,0,0]); %  glucose, water, oxygen

model = setParam(model,'obj',{'biomass_components'},1); %% Optimization of the biomass %%
model  = simplifyModel(model);

Sol = solveLP(model,'max'); 
printFluxes(model,Sol.x,true,10^-6,[],'%eqn\t%flux\n'); %% only exchange fluxes

idl =  find(ismember(string(model.subSystems),'Glycolysis / Gluconeogenesis'));
[num2cell(Sol.x(idl)), model.rxns(idl),constructEquations(model,model.rxns(idl)),string(model.subSystems(idl))]
Sol

disp('Recommended.. open excel for each models and do a manual curation ....')

return
                                %%%%%%%%%%%%%%%%%%%%% The End %%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
                                
                                                                   
% disp("Check tasks .... ");
% [taskReport, essentialRxns, taskStructure]=checkTasks(model_out,'/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Scripts/Metabolic_Task2.xlsx',1,0);          
% disp('done!')

%disp("Fill gaps .... ");
%[newConnected,cannotConnect,~,model_out, exitFlag] = fillGaps(REF{1},ref_model,false,true);


% exportModel(HTimmR,'URL',false); %% in RAVEN format
% exportToExcelFormat(HTimmR,'URL')

% % %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% % %%%%%% iMAT based contextualization %%%%%% Note: All fluxes were re-constrianted by this functions
% % %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % pp=1;
% %  disp("iMAT Contextualizing the liver model ? popularly used!");
% % 
% %    if pp
% %    tic
% %      [expressionRxns, parsedGPR] = mapExpressionToReactions(HTimmR,expression_normal); %% Generate reaction scores for the Healthy Obese/normal %%%
% %      options=struct();
% %      options.solver='iMAT';
% %      options.threshold_lb = 0;
% %      options.threshold_ub = 1.1;  %% Atleast 2 times expressed in the absolute scale
% %      options.expressionRxns = expressionRxns;
% %      
% %      HTimmR_normal_iMAT = createTissueSpecificModel(HTimmR, options);
% %      HTimmR_normal_iMAT.id = 'HMR_2.00_iMAT_normal';
% %      %save('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/Models/Contextualized/HTimmR_normal_iMAT','HTimmR_normal_iMAT');
% %    toc
% %       disp('done!')
% %    end
   
   