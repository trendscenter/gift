function fieldProp = icatb_extract_struct(structVar)

% input: structure array
% output: field names and values

if nargin > 1
    error('Only one argument needs to be passed');
end


if prod(size(struct)) > 1
    error('Input must be of size 1 by 1');
end

if ~isstruct(structVar)
    error('input must be of type structure');
end

% initialise output
fieldProp = cell(1, 2*length(structVar));

% extract the field names
fieldNames = fieldnames(structVar);

% loop over all fields
for ii = 1:length(fieldNames)
    fieldProp{2*ii - 1} = fieldNames{ii};
    fieldProp{2*ii} = getfield(structVar, fieldNames{ii});
end