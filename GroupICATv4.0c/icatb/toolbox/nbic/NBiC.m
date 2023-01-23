function bic = NBiC(x_size, y_size, tolerance, variable_Ids, reps, raw_data_matrix)
%%
% N-BiC is a biclustering method for simultaneously clustering two interacting variables
% It clusters a 2D matrix where each dimension represents a variable. Each variable has a set of examples 
% Variables could be subjects, neuro components
% It requires a sorting method to build the set of instancs for all given examples of the variable 
% If we have 'n' components
% After sorting, each component becomes a susbet of subjects 
% Prerequisite for the N-Bic method: Why sorting? 
% Since we use intersection for the merging or composing the biclusters,
% We need to represent variables from one dimension as a function (subset) of variables from other dimension 
% We sort the components as a subset of subejects based on a predefiend heuristic
% 'Search_BIC' is used to explore the biclusters for a given set of variables 
% The validator is integrated with the 'Search_BIC' subroutine
% A bicluster is a two-dimensional submatrix
% Input:
% x_size: minimum number fo rows (instances in x dimension) #components in a  biclusters
% y_size: minimum number fo cols (instances in y dimension) #subjects in a  biclusters
% tolerance: Allowed percentage of overlap between two biclusters 
% variable_Ids: Ids of sorted variables under consideration. If you want to run the analysis on 
% all the features, then it would just be 1:number of features. 
% But if we want to skip some features from the analysis, this option will provide us the flexibility   
% reps = no of expected repititions. The input ids would be randomly permutated 'reps' number of times  
% raw_data_matrix: the raw data
% Output:
% a list of biclusters 
%%
global minCmp;
minCmp = y_size;
global tol;
tol = tolerance;
global minSub;
minSub = x_size;
global BicId;
BicId = 1;
initBicList();
global BicList;
initListofBics();
global ListofBics;
%% Sorting or the features reduction
sorted_data = Feature_reduction(raw_data_matrix);
components = sorted_data;
%% Set of instances to bicluster 
% A set of Ids for variables under consideration. Each set is a vector of elements
% variable_Ids = 1:size(raw_data_matrix,2);
SyncComps = variable_Ids; % i.e., [1 5 7 13 14 16 17 28 30]; 
ComPlusSm = SyncComps;
permutations = perms(ComPlusSm);
%% Searching Biclusters for Intial run
simC = [];
simS = [];
if(length(ComPlusSm)>=minCmp)
Search_BIC(ComPlusSm,components,simC, simS,1);
ListofBics = BicList;
initBicList();
BicId = 1;
fprintf('Permuted Run 1 Completed! \n');
end
%end
%% Searching biclusters for secondary run: for other permutations depending on repititions expected  
% The validator would scrutinize before we add new biclusters from these runs
temp = randperm(length(permutations));
n_perm= temp(1:reps);
for perm = 1:reps 
if(length(n_perm(perm))>=minCmp)
fprintf('Permutation %u running... \n', perm);
Search_BIC(permutations(n_perm(perm),:),components,simC, simS,1); 
end

% Stability Checking 
itol =  30; % inter run overlap tolerance. We can use a different threhodls for both 'tol' and 'itol'  
len = length(ListofBics);
for b = 1: length(BicList)
    
for j = 1:len % Becuase We don't want to compare with new added Bics  
    % More than '100-itol' percent similar so should be redundant 
    overlapped_sub = length(intersect(BicList(b).subs,ListofBics(j).subs));
    Xpected_sub_overlap = ((100-itol)/100) * max(length(ListofBics(j).subs),length(BicList(b).subs));
    overlapped_cmp = length(intersect(BicList(b).comps,ListofBics(j).comps));
    Xpected_cmp_overlap = ((100-itol)/100) * max(length(ListofBics(j).comps),length(BicList(b).comps)); 
    if (overlapped_sub >= Xpected_sub_overlap && overlapped_cmp >= Xpected_cmp_overlap)
    ListofBics(j).freq  = ListofBics(j).freq +1;
    BicList(b).freq = -1;
    end
end
if(BicList(b).freq ~= -1 && ~isempty(BicList(b).subs))
ListofBics(length(ListofBics)+1) = BicList(b); 
end
end
BicId = 1;
initBicList();
end
bic = ListofBics;
fprintf('All the biclusters has been enlisted\n');
end
%%