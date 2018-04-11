function icatb_check_path
% Check path

try
    
    % Check startup file
    checkStartup;
    
catch
    
end

% Check required paths for the toolboxes
checkReqdPaths;

% run only on windows
icatb_openMatlabServer;


function checkStartup
% check gift start up file

if (exist('gift_startup.m', 'file') == 2) || (exist('gift_startup.m', 'file') == 6)
    fprintf( 'Executing gift_startup ...\n' );
    % execute gift start up file
    gift_startup;
    
else
    
    % Get the full file path of the gift file
    giftPath = which('gift.m');
    % Folder location of the gift
    giftPath = fileparts(giftPath);
    
    % all directories on path
    pathstr = path;
    
    % Get the directories on path
    allDirs = strread(pathstr, '%s', 'delimiter', pathsep);
    
    if ~isempty(allDirs)
        
        [indices] = regexp(allDirs, 'icatb$');
        indices = good_cells(indices);
        matchedDirs = allDirs(indices);
        
        if length(matchedDirs) > 1
            error('Fix MATLAB path such that it has only one version of GIFT at a time');
        elseif length(matchedDirs) == 1
            if strcmpi(giftPath, matchedDirs{1})
                return;
            else
                error('Fix MATLAB path such that it has only one version of GIFT at a time');
            end
        end
        
    end
    
end
% end for adding path

function ind = good_cells(mycell)

if ~iscell(mycell)
    mycell = {mycell};
end

ind = cellfun('isempty', mycell);

% Good cells
ind = (ind == 0);

% Convert ind to row vector
ind = ind(:)';

function checkReqdPaths
% Check required paths

% Get the full file path of the gift file
giftPath = which('gift.m');
% Folder location of the gift
giftPath = fileparts(giftPath);

reqdDirs = str2mat(giftPath, strcat(giftPath, filesep, str2mat('icatb_analysis_functions', 'icatb_batch_files', 'icatb_display_functions', ...
    'icatb_helpManual', 'icatb_helper_functions', 'icatb_io_data_functions', 'icatb_mex_files', 'icatb_talairach_scripts', 'icatb_spm8_files', 'icatb_parallel_files', ...
    'icatb_mancovan_files', ['toolbox', filesep, 'eegiftv1.0c'], ['toolbox', filesep, 'icasso122'], ['toolbox', filesep, 'mancovan'], ['toolbox', filesep, 'Graphical_Lasso'], ...
    ['toolbox', filesep, 'mi'], ['toolbox', filesep, 'noisecloud'], ['toolbox', filesep, 'export_fig'])));

pathstr = path;

allDirs = strread(pathstr, '%s', 'delimiter', pathsep);

% Add required paths
for nDir = 1:size(reqdDirs, 1)
    currentPath = deblank(reqdDirs(nDir, :));
    check = strmatch(currentPath, allDirs, 'exact');
    if isempty(check)
        addpath(genpath(currentPath), '-end');
    end
end
% End for adding required paths