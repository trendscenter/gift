function sorted_data = Feature_reduction(raw_data_matrix)
%% The given data matrix has a dimension of #samples x #features 
% We perform the feature reduction to represent each feature as a (reduced) subset of samples 
% Another rationale is to represent variables from one dimension as a function of variables from other dimension 
% Here, we sort the components (features) as a subset of subejects based on a predefiend heuristic. 
% We presented several approaches to perform the sorting in the referred paper.
% However, the sorting is more intuitively carried out based on the objectives of a study. 
% In another word, study dependent
%% The default method for feature reduction
% no_of_samples = size(raw_data_matrix,1);
no_of_features = size(raw_data_matrix,2);
sorted_data = cell(no_of_features, 1);
for i = 1:no_of_features
    comp = raw_data_matrix(:,i);
    pos_mean = mean(comp(comp>0));
    neg_mean = mean(comp(comp<0));
    poscomps = find(comp>=pos_mean);
    negcomps = find(comp<=neg_mean);
    sorted_data{i,1} = union(poscomps, negcomps); 
end
end