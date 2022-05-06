function [OverSubs] = ReporterMetaboliteEstimate(model,Genes,Pval1,FC1,baseOutURL,CondiStr,sigLevel)

ReporterFileOutURL= strcat(baseOutURL,CondiStr,'.txt'); %% construct
ReporterFileOutURL
EnrichmentFileOutURL= strcat(baseOutURL,'reporterM_Enriched','.xlsx'); %% construct
EnrichmentFileOutURL

% ** Works with the Cobra Structure as well! ** %
  disp(strcat('Performing reporter Metabolites for ',CondiStr, 'plz wait...'));
  
  repMets=reporterMetabolites(model,Genes,Pval1,true,ReporterFileOutURL,FC1);
     repMets_All = repMets(1);
     repMets_Up = repMets(2);
     repMets_Down = repMets(3);
     
  disp('done!');


  % OR analysis %
     
  disp(strcat('Performing overrepresentation analysis for ',CondiStr, 'plz wait...'));
  
     MetsID = [repMets_Up.mets(find(repMets_Up.metPValues < sigLevel));repMets_Down.mets(find(repMets_Down.metPValues < sigLevel))];
      [subs,iuu,Urxnsubs,RxnsSub,MetSub] = Subsystem(model,MetsID); %% [subs1';num2cell(iuu1)] Subsystem matched reactions
      %D=FEA(model,Urxnsubs1(:,1), char(string(Urxnsubs1(:,2))) ) 
      OR  = FEA(model,find(ismember(model.rxns,Urxnsubs(:,1))), 'subSystems');
      bix = find(str2double(string(OR(:,2))) > 0 & str2double(string(OR(:,2))) < sigLevel); 
      A = OR(bix,:);
      A = [A,num2cell((cell2mat(A(:,4))./cell2mat(A(:,5)))*100)];
      OR  = [[OR(1,:),'percentage']; A];
      xlwrite(EnrichmentFileOutURL,OR,CondiStr);
     MetsID =[];
     
     OverSubs = OR;
     disp('assign objects')
     disp('done!');
     disp(strcat('Results saved in .. ',ReporterFileOutURL));

end