function [input_data_file_patterns, result_bids] = icatb_parseBIDS(bids_info)
%% Parse BIDS
%

disp('Getting files information from BIDS ...');

modality_dir = 'func';
try
    modality_dir = bids_info.modality_dir;
catch
end

subjects = '';
try
    subjects = bids_info.subjects;
catch
end

sessions = '';
try
    sessions = bids_info.sessions;
catch
end

tasks = '';
try
    tasks = bids_info.tasks;
catch
end

result_bids = icatb_spm_BIDS(bids_info.root_dir);


%% Check subjects
if (~isempty(subjects))
    iab = ismember(lower(cellstr(char(result_bids.subjects.name))), lower(subjects));
    result_bids.subjects = result_bids.subjects(iab);
    %[iab, ia, ib] = intersect(lower(cellstr(char(result_bids.subjects.name))), lower(subjects));
    if (isempty(result_bids.subjects))
        error('Please check the subjects passed in the bids_info variable');
    end
    %result_bids.subjects = result_bids.subjects(ia);
end

%% Check sessions
sessionsIn = zeros(1, length(result_bids.subjects));
for nR = 1:length(sessionsIn)
    
    if (~isempty(sessions))
        ses_name =  result_bids.subjects(nR).session;
        chk = strmatch(lower(ses_name), sessions, 'exact');
        if (~isempty(chk))
            sessionsIn(nR) = 1;
        end
    else
        sessionsIn(nR) = 1;
    end
    
end


sessionsIn = find(sessionsIn == 1);
if (isempty(sessionsIn))
    error('Please check the sessions passed in the bids_info variable');
end

result_bids.subjects = result_bids.subjects(sessionsIn);

%% Filter based on tasks
if (strcmpi(modality_dir, 'func'))
    if (~isempty(tasks))
        
        subjectsIn = 0;
        for nR = 1:length(result_bids.subjects)
            
            all_tasks = cellstr(char(result_bids.subjects(nR).func.task));
            
            [iab, ia, ib] = intersect(lower(all_tasks), lower(tasks));
            
            subjectsIn = subjectsIn + length(ia);
            
            result_bids.subjects(nR).func = result_bids.subjects(nR).func(ia);
            
        end
        
        %end
    else
        subjectsIn = length(result_bids.subjects);
    end
end

if (subjectsIn == 0)
    error('No subjects found with the specified parameters in bids_info variable');
end


input_data_file_patterns = cell(length(result_bids.subjects), 1);

for nR = 1:length(result_bids.subjects)
    tmp = char(result_bids.subjects(nR).func.filename);
    tmp = icatb_fullFile('directory', fullfile(result_bids.subjects(nR).path, 'func'), 'files', tmp);
    input_data_file_patterns{nR} = cellstr(tmp);
end

input_data_file_patterns = cat(1, input_data_file_patterns{:});

chk = regexpi(input_data_file_patterns, '.*\.nii$|.*\.nii\.gz$|.*\.img$');
inds = icatb_good_cells(chk);
input_data_file_patterns = input_data_file_patterns(inds);

if (isempty(input_data_file_patterns))
    error('No fmri files found with the matching tasks');
end

disp('Done');
fprintf('\n');