function icatb_error(errorString, infoCell)
% icatb_error(errorString, infoString)
% errorString = string to be displayed via error function
% infoCell is a cell that contains information about the error to print

if(~exist('infoCell','var'))
    infoCell{1} ='';
end

for i = 1:size(infoCell,2)
    if(i==1)
        disp('-Error Information');
    end
    disp(infoCell{i});
end

error(errorString);