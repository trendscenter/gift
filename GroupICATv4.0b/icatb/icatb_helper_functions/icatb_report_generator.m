function icatb_report_generator(param_file, results)
%% Report generator for fmri and smri
%


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
    outDir = fullfile(fileparts(param_file), [sesInfo.userInput.prefix, '_gica_results']);
    opts.codeToEvaluate = 'icatb_gica_html_report(param_file, results);';
else
    resultsFile = 'icatb_sbm_html_report';
    outDir = fullfile(fileparts(param_file), [sesInfo.userInput.prefix, '_sbm_results']);
    opts.codeToEvaluate = 'icatb_sbm_html_report(param_file, results);';
end
opts.outputDir = outDir;
opts.showCode = false;
opts.useNewFigure = false;
opts.format = lower(formatName);
opts.createThumbnail = true;
if (strcmpi(opts.format, 'pdf'))
    opt.useNewFigure = false;
end
assignin('base', 'param_file', param_file);
assignin('base', 'results', results);

%publish('icatb_gica_html_report', 'outputDir', outDir, 'showCode', false, 'useNewFigure', false);
disp('Generating reults summary. Please wait ....');
drawnow;

publish(resultsFile, opts);

close all;

if (strcmpi(opts.format, 'html'))
    icatb_openHTMLHelpFile(fullfile(outDir, [resultsFile, '.html']));
else
    open(fullfile(outDir, [resultsFile, '.pdf']));
end

disp('Done');
