function Bic = Search_BIC(CMP_idx,components,simC,simS,start)

% Search biclusters by building different subsets of components
% It uses a modified depth-first-search (DFS) technique to create all possible subset of a given set
% And, traverse through the branches
% Use early abondoning the brances based on the imposed threshold
% The thresholds are usedefined
% User can control the size of the biclusters. Seperate for each dimension. Also, the ovrlap between them    
% Input:
% CMP_idx = An array of variables. It takes the Ids and the array seperately. 
% So that we can run the tweaks on the set of variables.   
% Ids: For example [1 3 5 8 9 6] 
% components = The cell array of instancs for a variable (sets).
% where each set is a group/vector of instances for a given variable
% simC = set of instances (for variable 1) in a BIC
% simS = set of instances (for variable 2) in a BIC
% start = A tracker of the component that being used to create the subset 
% by munna Dated: 4/Aug/2017

%% Write bicluster
 
global BicList;
global BicId;
global minSub;
global minCmp;
if(length(simS)>=minSub && length(simC)>=minCmp)    
if(BicId == 1) % Checking for first cluster
BicList(BicId).subs = simS;
BicList(BicId).comps = simC;
BicList(BicId).freq = 1;
%fprintf("biCluster# %u is done\n",BicId);      
BicId = BicId+1;
else
    %v = BIC_validation(simS,simC); 
    v = BIC_validation(simS,simC); % A different strategy for validation 
    if(v == 1)
     %fprintf("Adding biCluster# %u to the list....\n",BicId);      
     BicList(BicId).subs = simS;
     BicList(BicId).comps = simC;
     BicList(BicId).freq = 1;
     BicId = BicId+1;
    end
end
end
%% Generate Subsets of Components
for i = start:length(CMP_idx)
    if(components{CMP_idx(i),1}~= 0)
    if(isempty(simC))
        % Handle the base case!! Do something at first iteration of simCmp. 
        % The whole thing is just for initializing the comm 
        comm = [];
        comm = components{CMP_idx(i),1};
        %fprintf("ENTERED IN fOR LOOP\n",i);
        %simC(end+1)= CMP_idx(i);
        %comm = union(comm,components{simC(end),1});
        %simC(end) = [];
    else
        comm = components{simC(1),1}; % Taking the subjects of first component from the comps list 
        for j = 1:length(simC)
        comm = intersect(comm,components{simC(j),1});
        %fprintf("Secondary Iteration\n");
        end
    end
    %fprintf("Subsets Inside %u ....\n",simC);
    simS = intersect(comm,components{CMP_idx(i),1});
    if(length(simS)>=minSub)
        %fprintf("Enough Subjects\n");
        simC(end+1)= CMP_idx(i);
        %fprintf("Enough components %u\n",length(simC));
        Search_BIC(CMP_idx,components,simC,simS,i+1);
        simC(end) = [];
    end
    end
end
%% End of recursion
end