function mancovanInfo = icatb_mancovan_interactions(mancovanInfo, interactions)
%% Model Interactions

covariateTypes = cellstr(char(mancovanInfo.userInput.cov.type));
covariateNames = cellstr(char(mancovanInfo.userInput.cov.name));

group_ind = strmatch('categorical', lower(covariateTypes), 'exact');
cov_ind = strmatch('continuous', lower(covariateTypes), 'exact');
groups = covariateNames(group_ind);
covariates = covariateNames(cov_ind);

interactionsList = struct;
nS = 0;

if (length(group_ind) > 1)
    
    gglist = nchoosek((1:length(group_ind)), 2);
    
    for nG = 1:size(gglist, 1)
        % Don't show the interactions for multi-level interactions with
        % levels more than or equal to 3.
        levels = [length(unique(mancovanInfo.userInput.cov(gglist(nG, 1)).value)), length(unique(mancovanInfo.userInput.cov(gglist(nG, 2)).value))];
        if all(levels >= 3)
            disp(['Not Showing ', groups{gglist(nG, 1)}, ' X ', groups{gglist(nG, 2)}, ' interactions in the list as each has more than or equal to 3 levels']);
        else
            nS = nS + 1;
            interactionsList(nS).name = [groups{gglist(nG, 1)}, ' X ', groups{gglist(nG, 2)}];
            interactionsList(nS).type = 'group-group';
            interactionsList(nS).val = gglist(nG, :);
        end
        
    end
    
end

if (length(cov_ind) > 1)
    cclist = nchoosek((1:length(cov_ind)), 2);
    for nC = 1:size(cclist, 1)
        nS = nS + 1;
        interactionsList(nS).name = [covariates{cclist(nC, 1)}, ' X ', covariates{cclist(nC, 2)}];
        interactionsList(nS).type = 'covariate-covariate';
        interactionsList(nS).val = cclist(nC, :);
    end
end

for i = 1:length(group_ind)
    for j = 1:length(cov_ind)
        nS = nS + 1;
        interactionsList(nS).name = [groups{i}, ' X ', covariates{j}];
        interactionsList(nS).type = 'group-covariate';
        interactionsList(nS).val = [i, j];
    end
end

sel = [];
if (isfield(interactionsList, 'name'))
    allInterNames = cellstr(char(interactionsList.name));
    if (~exist('interactions', 'var'))
        sel = icatb_listdlg('PromptString', 'Select Model Interactions ...', 'title_fig', 'Select Model Interactions ...', 'SelectionMode', 'multiple', 'ListString', allInterNames, 'movegui', 'center', 'windowStyle', 'modal', 'help', struct('title', 'Model Interactions', 'str', ...
            'Select model interactions which could be categorical-categorical, categorical-continuous, continuous-continuous. If you don''t want to include interactions, click on cancel button'));
    else
        for nInter = 1:size(interactions, 1)
            cName1 = [covariateNames{interactions(nInter, 1)}, ' X ', covariateNames{interactions(nInter, 2)}];
            chkIndex = strmatch(cName1, allInterNames, 'exact');
            if (isempty(chkIndex))
                cName2 = [covariateNames{interactions(nInter, 2)}, ' X ', covariateNames{interactions(nInter, 1)}];
                chkIndex = strmatch(cName2, allInterNames, 'exact');
            end
            if (~isempty(chkIndex))
                sel = [sel, chkIndex];
            end
        end
    end
end

drawnow;

gglist = [];
cclist = [];
gclist = [];
selLists = {};

if (~isempty(sel))
    
    selInteractionsList = interactionsList(sel);
    
    disp('Selected model interactions are: ');
    fprintf('\n');
    disp(char(selInteractionsList.name));
    fprintf('\n');
    
    gg_inds = strmatch('group-group', lower(cellstr(char(selInteractionsList.type))), 'exact');
    cc_inds = strmatch('covariate-covariate', lower(cellstr(char(selInteractionsList.type))), 'exact');
    gc_inds = strmatch('group-covariate', lower(cellstr(char(selInteractionsList.type))), 'exact');
    
    gglist = cat(1, selInteractionsList(gg_inds).val);
    cclist = cat(1, selInteractionsList(cc_inds).val);
    gclist = cat(1, selInteractionsList(gc_inds).val);
    
    if (~isempty(gg_inds))
        selLists{end + 1} = 'group-group';
    end
    
    if (~isempty(cc_inds))
        selLists{end + 1} = 'covariate-covariate';
    end
    
    if (~isempty(gc_inds))
        selLists{end + 1} = 'group-covariate';
    end
    
    mancovanInfo.userInput.modelInteractions.list = cellstr(char(selInteractionsList.name));
    
end

mancovanInfo.userInput.modelInteractions.gglist = gglist;
mancovanInfo.userInput.modelInteractions.cclist = cclist;
mancovanInfo.userInput.modelInteractions.gclist = gclist;
mancovanInfo.userInput.modelInteractions.types = selLists;