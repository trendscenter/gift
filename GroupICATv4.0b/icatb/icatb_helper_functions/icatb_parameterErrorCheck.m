function icatb_parameterErrorCheck(sesInfo, optional)
%  icatb_parameterErrorCheck(sesInfo)
% function to make sure parameters enter by user make sense

if ~exist('optional', 'var')
    optional = 'print';
end


if strcmpi(optional, 'print')
    
    disp('');
    disp('Checking to make sure parameters are correct...');
    
end

[modalityType, dataTitle] = icatb_get_modality;

%make sure all file sets are the same size
numOfFileSets = length(sesInfo.inputFiles);

%numOfScans = size(sesInfo.inputFiles(1).name,1);
% get the number of time points
[numOfScans] = icatb_get_countTimePoints(sesInfo.inputFiles(1));

% making it to work for different time points
% for fileSet = 2:numOfFileSets
%     numOfScans2 = size(sesInfo.inputFiles(fileSet).name,1);
%     if( numOfScans2 ~= numOfScans)
%         infoCell{1} = sesInfo.inputFiles(fileSet).name;
%         infoCell{2} = [' The listed files above contains ',num2str(numOfScans2), ' files while the first file set contains ',num2str(numOfScans)];
%         icatb_error('SELECTED FILES DON''T MATCH:',infoCell);
%     end
% end

%make sure there are voxels in mask
if strcmpi(optional, 'print')
    disp('  Checking mask ');
end

if(isempty(sesInfo.mask_ind))
    %icatb_error(['BRAIN MASK CONTAINS NO VOXELS']);
    if ~strcmpi(modalityType, 'eeg') 
        error('Brain mask contains no voxels');
    else
        error('EEG mask contains zero indices');
    end
end

%make sure number of pc's make sense.
if strcmpi(optional, 'print')
    disp('  Checking principal component parameters');
end

% if(numOfScans < sesInfo.reduction(1).numOfPCAfterReduction)
%     infoCell{1} = 'Reduction Step 1: ';
%     infoCell{2} = ['Maximum number allowed for PC1 is ', num2str(numOfScans)];
%     infoCell{3} = ['Number you selected for PC1 is ', num2str(sesInfo.reduction(1).numOfPCAfterReduction)];
%     errorString = 'More principal components than timepoints.';
%     %error(['More principal components than timepoints']);
%     icatb_error(errorString, infoCell);
% end

for i=1:sesInfo.numReductionSteps
    if(sesInfo.reduction(i).numOfPCInEachGroupAfterCAT(1) < sesInfo.reduction(i).numOfPCAfterReduction)
        infoCell{1} = ['Reduction Step ', num2str(i), ': '];
        infoCell{2} = ['Maximum number allowed for PC', num2str(i), ' is ', ...
                num2str(sesInfo.reduction(i).numOfPCInEachGroupAfterCAT(1))];
        infoCell{3} = ['Number you selected for PC', num2str(i), ' is ', ...
                num2str(sesInfo.reduction(i).numOfPCAfterReduction(1))];
        errorString = 'More principal components than timepoints.';
        icatb_error(errorString, infoCell);
    end
end


%make sure header files match


if strcmpi(optional, 'print')
    disp('Done with parameter error check ');
    disp('');
end