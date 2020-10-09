function R2 = icatb_model_order_prediction(obs_file_name, ref_file_name, varargin)
%% Model order prediction
% Multiple regression framework is used to determine variances contributed
% by higher model order components to the lower model order components.
%
%
% Inputs:
% 1. obs_file_name - Lower model order components in a nifti file
% 2. ref_file_name - Higher model order components in a nifti file
% 3. varargin - Arguments passsed in pairs
%   a. threshold - Variance threshold.
%

icatb_defaults;
global FONT_COLOR;

%% Select lower model order components
if (~exist('obs_file_name', 'var'))
    obs_file_name = icatb_selectEntry('typeSelection', 'multiple', 'typeEntity', 'file', 'filter', '*agg*comp*.img;*agg*comp*.nii', 'title', 'Select lower model order components', ...
        'filetype', 'image');
    drawnow;
end

if (isempty(obs_file_name))
    error('Lower model order components are not selected');
end


%% Select higher model order components
if (~exist('ref_file_name', 'var'))
    ref_file_name = icatb_selectEntry('typeSelection', 'multiple', 'typeEntity', 'file', 'filter', '*agg*comp*.img;*agg*comp*.nii', 'title', 'Select higher model order components', ...
        'filetype', 'image');
    drawnow;
end

if (isempty(ref_file_name))
    error('Higher model order components are not selected');
end

threshold = 10;
displayResults = 0;


if (~isempty(varargin))
    
    for i = 1:2:length(varargin)
        if (strcmpi(varargin{i}, 'threshold'))
            threshold = varargin{i + 1};
        elseif (strcmpi(varargin{i}, 'display'))
            displayResults = varargin{i + 1};
        end
    end
    
else
    
    dlg_title = 'Select options for lower vs higher model order';
    
    numParameters = 1;
    
    inputText(numParameters).promptString = 'Enter variance threshold in %';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(threshold);
    inputText(numParameters).dataType = 'numeric';
    inputText(numParameters).tag = 'threshold';
    inputText(numParameters).enable = 'on';
    
    answers = icatb_inputDialog('inputtext', inputText, 'Title', dlg_title, 'handle_visibility', 'on');
    
    if (~isempty(answers))
        threshold = answers{1};
    end
    
end


obs_file_name = icatb_rename_4d_file(obs_file_name);
ref_file_name = icatb_rename_4d_file(ref_file_name);

nfiles = size(obs_file_name, 1);

%% Get R^2 summary (R2 - No. of comps in higher model order x no. of comps in lower model order)
disp('Computing variance captured by higher model order using lower model order components as observations ...');
strs = repmat({''}, nfiles, 1);
for i = 1:nfiles
    
    y = icatb_read_data(deblank(obs_file_name(i, :)));
    maskN = find(abs(y) > eps);
    y = icatb_remove_mean(y(maskN));
    ref = icatb_remove_mean(icatb_read_data(ref_file_name, [], maskN));
    betas = pinv(ref)*y;
    denom = sum(abs(eig(icatb_cov(y))));
    r2 = zeros(1, length(betas));
    for nR = 1:length(betas)
        res = y - betas(nR)*ref(:, nR);
        num = sum(abs(eig(icatb_cov(res))));
        r2(nR) = 100*(1 - num/denom);
    end
    
    if (i == 1)
        R2 = zeros(length(r2), nfiles);
    end
    
    R2(:, i) = r2;
    
    nums = find(abs(r2) > threshold);
    vals = r2(nums);
    str = '';
    for nval = 1:length(vals)
        str = [str, num2str(nums(nval)), '(', num2str(vals(nval), '%0.3f'), '%),'];
    end
    str(end) = '';
    strs{i} = str;
end

clear str;

if (displayResults)
    %% Write variance info of dimensions higher model order x lower model order
    [pathstr, fN, extn] = fileparts(ref_file_name(1, :));
    outFileName = fullfile(pathstr, [fN, '_lm_vs_hm.txt']);
    
    numPara = 1;
    varStruct(numPara).tag = 'Lower Model Order Component #';
    varStruct(numPara).value = (1:nfiles)';
    
    numPara = numPara + 1;
    varStruct(numPara).tag = 'Higher Model Order Component # (variance in %)';
    varStruct(numPara).value = char(strs);
    
    icatb_printToFile(outFileName, varStruct, 'Lower vs Higher Model order Summary', 'column_wise');
    
    
    %dlmwrite(outFileName, R2, 'precision', '%0.2f');
    disp(['Variance info is saved in file ', outFileName]);
    %disp('R^2 is of dimensions higher model order comps x lower model order comps');
    fprintf('\n');
    
    R2(abs(R2) < threshold) = 0;
    
    
    %% Display images
    xticks = cellstr(num2str((1:size(R2, 2))'));
    yticks = cellstr(num2str((1:size(R2, 1))'));
    [aH, cH] = icatb_plot_matrix(R2', xticks, yticks, 'cmap', jet(64));
    sh = get(aH, 'currentAxes');
    xlabel('Lower Model Order', 'parent', sh);
    ylabel('Higher Model Order', 'parent', sh);
    ylabel('Variance (%)', 'parent', cH);
    set(sh, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);
    set(cH, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);
end



