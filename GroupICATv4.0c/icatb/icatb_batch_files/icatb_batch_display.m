function icatb_batch_display(inputFile)
% Batch script for displaying component images
% Input file must be in fullfile path

[pathstr, fName, extn] = fileparts(inputFile);

if isempty(pathstr)
    tempF = inputFile;
    inputFile = which(inputFile);
    if isempty(inputFile)
        error('Error:InpFile', 'File %s doesn''t exist', tempF);
    end
end

[parameters] = icatb_readParameters_display(inputFile);

if strcmpi(parameters.displayType, 'orthogonal viewer')
    icatb_orthoViewer([], parameters );
elseif strcmpi(parameters.displayType, 'composite viewer')
    icatb_compositeViewer([], parameters);
else
    icatb_componentExplore([], parameters);
end
