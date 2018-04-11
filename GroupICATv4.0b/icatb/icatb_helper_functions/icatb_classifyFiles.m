function [newPrefix, indices_prefix, originalPrefix] = icatb_classifyFiles(files)
% function to get the prefixes of the corresponding files
% and shows the number of files to the user

[prefix, indices_prefix] = icatb_spm_get('fileSummary', files);
originalPrefix = prefix;

prefixStruct = repmat(struct('prefix', ''), 1, size(prefix, 1));

% Show the number of entries during the prefix mode
for ii = 1:size(prefix, 1)
    temp = deblank(prefix(ii, :));
    findNumberFiles = find(indices_prefix(ii, :) ~= 0);
    temp = ['[1:', num2str(length(findNumberFiles)), '] ', temp];
    prefixStruct(ii).prefix = temp;    
end
clear prefix;
newPrefix = str2mat(prefixStruct.prefix); % Updated prefix