function icatb_report_generator(param_file, results)
%% Report generator for fmri and smri
%

sesDir = fileparts(param_file);

if (isempty(sesDir))
    sesDir = pwd;
end

load (param_file);

modalityType = 'fmri';

try
    modalityType = sesInfo.modality;
catch
end

appName = 'group_ica_modality';

setappdata(0, appName, modalityType);

try
    formatName = results.formatName;
catch
end

if (~exist('formatName', 'var'))
    formatName = 'html';
end

if (isnumeric(formatName))
    if (formatName == 1)
        formatName = 'html';
    else
        formatName = 'pdf';
    end
end

results.formatName = formatName;

if (strcmpi(modalityType, 'fmri'))
    resultsFile = 'icatb_gica_html_report';
    outDir = fullfile(sesDir, [sesInfo.userInput.prefix, '_gica_results']);
    opts.codeToEvaluate = 'icatb_gica_html_report(param_file, results);';
else
    resultsFile = 'icatb_sbm_html_report';
    outDir = fullfile(sesDir, [sesInfo.userInput.prefix, '_sbm_results']);
    opts.codeToEvaluate = 'icatb_sbm_html_report(param_file, results);';
end

results.outputDir = outDir;
opts.outputDir = outDir;
opts.showCode = false;
opts.useNewFigure = false;
opts.format = lower(formatName);
opts.createThumbnail = true;
if (strcmpi(opts.format, 'pdf'))
    opts.useNewFigure = false;
end


disp('Generating reults summary. Please wait ....');
drawnow;

isDeployedVer = isdeployed;

if (~isDeployedVer)
    
    assignin('base', 'param_file', param_file);
    assignin('base', 'results', results);
    
    %publish('icatb_gica_html_report', 'outputDir', outDir, 'showCode', false, 'useNewFigure', false);
    % disp('Generating reults summary. Please wait ....');
    % drawnow;
    
    publish(resultsFile, opts);
    
    close all;
    
else
    results.save_info = 1;
    if (strcmpi(modalityType, 'fmri'))
        icatb_gica_html_report(param_file, results);
    elseif (strcmpi(modalityType, 'smri'))
        icatb_sbm_html_report(param_file, results);
    end
    
end


if (strcmpi(modalityType, 'fmri'))
    try
        network_opts = results.network_summary_opts;
        if (~isempty(network_opts))
            compFileNaming = sesInfo.icaOutputFiles(1).ses(1).name;
            compFiles = icatb_rename_4d_file(icatb_fullFile('directory', sesInfo.outputDir, 'files', compFileNaming));
            network_opts.file_names = compFiles;
            postProcessFile = fullfile(sesInfo.outputDir, [sesInfo.userInput.prefix, '_postprocess_results.mat']);
            load(postProcessFile, 'fnc_corrs_all');
            if (sesInfo.numOfSess > 1)
                fnc_corrs_all = reshape(mean(fnc_corrs_all, 2), sesInfo.numOfSub, sesInfo.numComp, sesInfo.numComp);
            else
                fnc_corrs_all = reshape(squeeze(fnc_corrs_all), sesInfo.numOfSub, sesInfo.numComp, sesInfo.numComp);
            end
            if (sesInfo.numOfSub > 1)
                fnc_corrs_all = squeeze(mean(fnc_corrs_all));
            else
                fnc_corrs_all = squeeze(fnc_corrs_all);
            end
            fnc_corrs_all = icatb_z_to_r(fnc_corrs_all);
            network_opts.fnc_matrix_file = fnc_corrs_all;
            network_opts.save_info = 1;
            if (~isfield(network_opts, 'outputDir'))
                network_opts.outputDir = fullfile(sesInfo.outputDir, 'network_summary');
            end
            if (~isfield(network_opts, 'prefix'))
                network_opts.prefix = [sesInfo.userInput.prefix, '_network_summary'];
            end
            icatb_network_summary(network_opts);
        end
    catch
    end
    
end


if (~isDeployedVer)
    if (strcmpi(opts.format, 'html'))
        icatb_openHTMLHelpFile(fullfile(outDir, [resultsFile, '.html']));
    else
        open(fullfile(outDir, [resultsFile, '.pdf']));
    end
end

disp('Done');
