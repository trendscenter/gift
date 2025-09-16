function structVar = icatb_includeFields_struct(structVar, fieldsReq)
% includes only the required fields in a structure

% all the fields
allFields = fieldnames(structVar);
% initialise vars
kk = 0; fields_to_remove = {};
% loop over all fields
for ii = 1:length(allFields)
    matchedIndex = strmatch(allFields{ii}, fieldsReq, 'exact');
    % collect fields to be removed
    if isempty(matchedIndex)
        kk = kk + 1;
        fields_to_remove{kk} = allFields{ii};
    end
end

% remove the corresponding fields in the structure
if ~isempty(fields_to_remove)
    structVar = rmfield(structVar, fields_to_remove);
end