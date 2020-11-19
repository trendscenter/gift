function icatb_copy_ica_results_to_bids(sesInfo)
%% Copy ica results to bids derivatives
% IC maps - *ic_maps.nii
% Timecourses - *timecourses.nii
%

icatb_defaults;
global COMPONENT_NAMING;
global TIMECOURSE_NAMING;

bids_dir = sesInfo.userInput.bids_info.root_dir;
derivatives_dir = '';
try
    derivatives_dir = sesInfo.userInput.bids_info.derivatives_dir;
catch
end

if (isempty(derivatives_dir))
    derivatives_dir = fullfile(bids_dir, 'derivatives', 'gift');
end

if (exist(derivatives_dir, 'dir') ~= 7)
    mkdir(derivatives_dir);
end

ica_results_dir = sesInfo.outputDir;
subjectICAFiles = icatb_parseOutputFiles('icaOutputFiles', sesInfo.icaOutputFiles, 'numOfSub', sesInfo.numOfSub, 'numOfSess', sesInfo.numOfSess);

for nP = 1:length(subjectICAFiles)
    
    tmp_fname = deblank(sesInfo.inputFiles(nP).name(1, :));
    [inp_pth, inp_file, extn] = fileparts (tmp_fname);
    
    disp(['Writing subject ica results of dataset ', fullfile(inp_pth, inp_file), ' to bids derivatives']);
    
    %% Get components
    P = subjectICAFiles(nP).ses(1).name;
    P = icatb_parseExtn(deblank(P(1, :)));
    
    % image extension
    [pathstr_comp, file_comp, imExtn] = fileparts(P);
    lastUnderScore = icatb_findstr(deblank(P(1, :)),'_');
    lastUnderScore = lastUnderScore(end);
    component_name = deblank(P(1,1:lastUnderScore));
    timecourse_name = strrep(component_name, COMPONENT_NAMING, TIMECOURSE_NAMING);
    timecourse_file = [timecourse_name, imExtn];
    
    components_file = fullfile(ica_results_dir, P);
    timecourse_file = fullfile(ica_results_dir, timecourse_file);
    
    suffix = strrep(inp_pth, bids_dir, '');
    bids_results_dir = fullfile(derivatives_dir, suffix);
    
    if (exist(bids_results_dir, 'dir') ~= 7)
        mkdir(bids_results_dir);
    end
    
    inp_file = strrep(inp_file, '.nii', '');
    
    new_components_file = fullfile(bids_results_dir, [inp_file, '_ic_maps.nii']);
    new_timecourse_file = fullfile(bids_results_dir, [inp_file, '_timecourses.nii']);
    
    if (~strcmpi(imExtn, '.nii'))
        %% Save component maps to bids derivatives
        compV = icatb_spm_vol(components_file);
        compData = icatb_spm_read_vols(compV);
        icatb_write_nifti_data(new_components_file, compV, compData);
        
        %% Save timecourse to bids derivatives
        timecourse = icatb_loadData(timecourse_file);
        timeV = compV(1);
        timeV.dim = [size(timecourse), 1];
        timeV.fname = new_timecourse_file;
        icatb_write_vol(timeV, timecourse);
    else
        copyfile(components_file, new_components_file);
        copyfile(timecourse_file, new_timecourse_file);
        
    end
    
end

fprintf('Done\n');