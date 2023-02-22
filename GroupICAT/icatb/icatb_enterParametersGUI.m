function icatb_enterParametersGUI(sub_file, inputFile)
%% Set up analysis


if exist('sub_file', 'var')
    icatb_setup_analysis(sub_file);
    return;
end


if exist('inputFile', 'var')
    icatb_setup_analysis(inputFile);
    return;
end


icatb_setup_analysis;