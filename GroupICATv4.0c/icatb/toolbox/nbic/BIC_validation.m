function validity = BIC_validation(simSub,simCmp)
% Applied F1 measure/Dice index between the BICs to control the replication
% The validation scrutinizes for the overlap between already reported and newly genrated biclusters
% Across a single permutation 
% If we wanna use F1 index then, for 20% overlap, we can allow <=0.37 as the threshold
global tol;
global minSub;
global minCmp;
global BicId;
global BicList;
ovrlookedSub = (tol*minSub)/100; % That portion we wanna allow between two BICs 
ovrlookedCmp = (tol*minCmp)/100;
if(ovrlookedCmp<1)
  ovrlookedCmp = 1;
end
%validity = 1;
%% Check validity
Nj = length(simSub)*length(simCmp); % Size of the current BIC
maxF1 = -1; % Intialization of maximum F1 score
simBICidx = [];
duplicate = 0;
for cl =1:BicId-1 % For the first one it will not enter to this for loop
    
    % F measure for current BIC to the others already in the list
    Sjk = length(intersect(simSub,BicList(cl).subs));
    Cjk = length(intersect(simCmp,BicList(cl).comps));
    Nk = length(BicList(cl).subs)* length(BicList(cl).comps);
    
    % ---------------------------------------------------------------------------------------
    if(Sjk>=minSub && Cjk<=ovrlookedCmp) % If subjects are almost similar but no comps similar or <=1
        BicList(cl).subs = intersect(BicList(cl).subs,simSub);
        BicList(cl).comps = union(BicList(cl).comps,simCmp);
        duplicate = 1; % Flag for duplicate
        BicList(cl).freq = BicList(cl).freq+1;
        fprintf('similar voxels\n')
        break;
    end
    if(Cjk>=minCmp && Sjk<=ovrlookedSub) % If Components are almost similar but no Subjects similar or <=1
        BicList(cl).subs = union(BicList(cl).subs,simSub);
        BicList(cl).comps = intersect(BicList(cl).comps,simCmp); 
        %BicList(cl).comps = union(BicList(cl).comps,simCmp);% Changed very recently for more traits  
        duplicate = 1; % Flag for duplicate
        BicList(cl).freq = BicList(cl).freq+1;
        fprintf('similar comps\n')
        break;
    end
    F1 = (2*Sjk*Cjk)/(Nj+Nk);
    if(F1>maxF1)
        maxF1 = F1;
        simBICidx = cl;
    end   
end
% After comparing the current with all BICs from the BicList we will definitely get a max F1
% Now if it crosses the threshold merge the BIC with older one and invalidate for not being duplicated  
% otherwise, validate to add it onto the BicList 
if(duplicate==1)
    validity = 0;    
else
if(maxF1>=0.28) 
% Increase the threshold to grab more biclusters [Goto else block] 
   
  % 1.--------------------- Merge by union ------------------------------------------- 
    %BicList(simBICidx).subs = union(BicList(simBICidx).subs,simSub);
    %BicList(simBICidx).comps = union(BicList(simBICidx).comps,simCmp);
    
  % 2.--------------------- Merge by intersection -------------------------------
  interSub = intersect(BicList(simBICidx).subs,simSub);
  interCmp = intersect(BicList(simBICidx).comps,simCmp);
  if(length(interSub)>=minSub && length(interCmp)>=minCmp)
  BicList(simBICidx).subs = interSub;
  BicList(simBICidx).comps = interCmp;
  fprintf('Merged in %d!\n',simBICidx);
  end
  %--------------------------------------------------------------- 
  % No new BIC, so, Don't increase the BicId
  %fprintf('Frequency increased at %d!\n',simBICidx);
  BicList(simBICidx).freq = BicList(simBICidx).freq+1;
  validity = 0;
else
  validity = 1;
end
end
end
