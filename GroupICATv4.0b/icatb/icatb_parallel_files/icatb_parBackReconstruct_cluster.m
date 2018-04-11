function varargout = icatb_parBackReconstruct_cluster(sesInfo, dataSetsToRun, verbose)
%% Back-reconstruct components in parallel mode for spatial temporal regression, moo-icar and temporal ica only
%

if (ischar(sesInfo))
    load(sesInfo);
end

if (~exist('sesInfo', 'var'))
    error('ICA Parameter file is not specified');
end

if (~exist('verbose', 'var'))
    verbose = 0;
end


inputFiles = sesInfo.inputFiles;
mask_ind = sesInfo.mask_ind;
preproc_type = sesInfo.preproc_type;
backReconType = sesInfo.backReconType;
back_reconstruction_mat_file = sesInfo.back_reconstruction_mat_file;
outputDir = sesInfo.outputDir;
diffTimePoints = sesInfo.diffTimePoints;
numOfSub = sesInfo.numOfSub;
numOfSess = sesInfo.numOfSess;

try
    modalityType = sesInfo.modality;
    setappdata(0, 'group_ica_modality', modalityType);
catch
end

[modalityType, dataTitle, compSetFields] = icatb_get_modality;

useTemporalICA = 0;
if (strcmpi(modalityType, 'fmri'))
    try
        useTemporalICA = strcmpi(sesInfo.group_ica_type, 'temporal');
    catch
    end
end

icain =[sesInfo.ica_mat_file, '.mat'];

if (~useTemporalICA)
    load(fullfile(outputDir, icain), 'icasig');
    numComp = size(icasig, 1);
    if (strcmpi(backReconType, 'spatial-temporal regression'))
        icasig = icasig';
    end
else
    load(fullfile(outputDir, icain), 'temporal_icasig');
    icasig = temporal_icasig;
    clear temporal_icasig;
    numComp = size(icasig, 1);
    icasig = icatb_remove_mean(icasig');
end

meanMap = zeros(numComp, length(mask_ind));

%% Loop over files
parfor nDataSet = 1:length(dataSetsToRun)
    
    tmpDataSetsToRun = dataSetsToRun;
    nSet = tmpDataSetsToRun(nDataSet);
    
    tmpFiles = inputFiles;
    
    tmpTp = diffTimePoints;
    
    tmpCompSetFields = compSetFields;
    
    if (verbose)
        disp(['Back reconstructing  set ', num2str(nSet)]);
    end
    
    % Read data
    data = icatb_read_data(char(tmpFiles(nSet).name), [], mask_ind);
    
    % Call pre-processing function
    data = icatb_preproc_data(data, preproc_type, verbose);
    
    drawnow;
    
    if (~useTemporalICA)
        
        if (strcmpi(backReconType, 'spatial-temporal regression'))
            %% STR
            if (numOfSub*numOfSess ~= 1)
                
                % Dual regression
                [tc, spatial_maps] = icatb_dual_regress(data, icasig);
                
            else
                
                %% For single subject ICA temporal regression is avoided
                tc = pinv(icatb_remove_mean(icasig))*icatb_remove_mean(data);
                spatial_maps = icasig';
                
            end
            
        else
            %% MOO-ICAR
            [spatial_maps, tc] = icatb_gigicar(data', icasig);
            
            if (numOfSub*numOfSess == 1)
                tc = tc';
            end
            
        end
        
    else
        %% Temporal ica
        data = data';
        if (nSet == 1)
            startT = 1;
        else
            startT = sum(tmpTp(1:nSet-1)) + 1;
        end
        endT = sum(tmpTp(1:nSet));
        
        tc = icasig(startT:endT, :);
        spatial_maps = pinv(tc)*data;
        
        if (numOfSub*numOfSess == 1)
            tc = tc';
        end
    end
    
    
    compSet = struct(tmpCompSetFields{1}, spatial_maps, tmpCompSetFields{2}, tc);
    
    meanMap = meanMap + spatial_maps;
    
    
    if (verbose)
        disp(['-done back reconstructing  set ', num2str(nSet)]);
    end
    
    drawnow;
    
    % Save Results
    subFile = [back_reconstruction_mat_file, num2str(nSet), '.mat'];
    
    if (verbose)
        txtMsg = ['-saving back reconstructed ica data for set ', num2str(nSet), ' -> ',subFile];
        disp(txtMsg);
        drawnow;
        icatb_parSave(fullfile(outputDir, subFile), {compSet}, {'compSet'});
    end
    
end
% End of loop over files


meanMap = meanMap./length(dataSetsToRun);

if (nargout == 1)
    varargout{1} = meanMap;
else
    varargout = {};
end