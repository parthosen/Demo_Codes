function [xx] = simulate_PS_Obj(model,constrainON)

%[selExc, selUpt] = findExcRxns(model); % Cobra

if(constrainON)
[rxns_uptake, uptIDX] = getExchangeRxns(model,'in'); %% -1000 - 1000
[rxns_prod, prodIDX] = getExchangeRxns(model,'out'); %% 0 - 1000
[rxns_both, bothIDX] = getExchangeRxns(model,'both'); %% -1000 - 1000 %% All exchange reactions
 model = setParam(model,'eq',bothIDX,0); %% Cannot uptake but can secrete
 model = setParam(model,'ub',prodIDX,10000); %% Cannot uptake but can secrete
end

model = setParam(model,'ub', {'HMR_9034','HMR_9047','HMR_9048','HMR_9058'},[1000,1000,1000,1000]); % glucose, water, oxygen
model = setParam(model,'lb', {'HMR_9034','HMR_9047','HMR_9048','HMR_9058'},[-1000,-1000,0,0]); %  glucose, water, oxygen

model = setParam(model,'obj',{'biomass_components'},1); %% Optimization of the biomass %%

model = simplifyModel(model);

Sol = optimizeCbModel(model);
xx = Sol;

end