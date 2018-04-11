function [outputFiles] = icatb_getOutputFileNames(prefix, numComp, numSub, numSess, writeInDirs)
% returns the output file naming
%
% Input:
% 1. prefix - output prefix
% 2. numComp - number of components
% 3. numSub - number of subjects
% 4. numSess - number of sessions
% 5. flagTimePoints - flag for time points ('same' or 'different')
%
% Output:
%
% outputFiles - output file naming


% load defaults
icatb_defaults;
global MEAN_AN3_FILE;
global TMAP_AN3_FILE;
global STD_AN3_FILE;
global MEAN_ALL_AN3_FILE;
global SUBJECT_ICA_AN3_FILE;
global SESSION_POSTFIX;
global COMPONENT_NAMING;
global TIMECOURSE_NAMING;

global GROUP_ICA_INDEX;
global MEAN_INDEX;
global TMAP_INDEX;
global STD_INDEX;
global SUBJECT_ICA_INDEX;
global MEAN_ALL_INDEX;
% functional data filter (also the variable needed to write the images)
global FUNCTIONAL_DATA_FILTER;


[modalityType, dataTitle] = icatb_get_modality;

if ~strcmpi(modalityType, 'eeg')
    % image extension
    [pp, bb, fileExtn] = fileparts(FUNCTIONAL_DATA_FILTER);
    % end for image extension

else
    fileExtn = '.mat';
end

% assign same time points by default
if ~exist('writeInDirs', 'var')
    writeInDirs = 0;
end

subNaming = SUBJECT_ICA_AN3_FILE;
meanAllNaming = MEAN_ALL_AN3_FILE;
meanNaming = MEAN_AN3_FILE;
stdNaming = STD_AN3_FILE;
tmapNaming = TMAP_AN3_FILE;

if (writeInDirs == 1)
    subNaming = fullfile('_scaling_components_files', regexprep(subNaming, '^(_)', ''));
    meanNaming = fullfile('_group_stats_files', regexprep(meanNaming, '^(_)', ''));
    stdNaming = fullfile('_group_stats_files', regexprep(stdNaming, '^(_)', ''));
    tmapNaming = fullfile('_group_stats_files', regexprep(tmapNaming, '^(_)', ''));
    meanAllNaming = fullfile('_group_stats_files', regexprep(meanAllNaming, '^(_)', ''));
end


% following are the cases:
% 1. one subject and one session (no mean is calculated)
% 2. one subject and multiple sessions (mean is calculated over all
% sessions)
% 3. multiple subjects and multiple sessions (stats are calculated)
% 4. for different time points mean is not calculated


if (~strcmpi(modalityType, 'smri'))

    if numSub*numSess == 1
        % for single subject single session

        % check image extensions
        if strcmpi(fileExtn, '.img')
            % analyze data

            % loop over components
            for kk = 1:numComp
                % component index
                fileIndex = icatb_returnFileIndex(kk);
                % single subject single session component naming
                outputFiles(1).ses(1).name(kk, :) = ...
                    [prefix, subNaming, '01', COMPONENT_NAMING, SESSION_POSTFIX, '1', ...
                    '_', fileIndex, fileExtn];

            end
            % end loop over components

        elseif strcmpi(fileExtn, '.nii') || strcmpi(fileExtn, '.mat')
            % nifti data

            % store all the components in one file
            outputFiles(1).ses(1).name = ...
                [prefix, subNaming, '01', COMPONENT_NAMING, SESSION_POSTFIX, '1', '_', fileExtn];
        else
            error('Unknown image extensions');
        end
        % end for checking image extensions

    elseif numSub == 1 && numSess > 1
        % one subject multiple sessions

        % check image extensions
        if strcmpi(fileExtn, '.img')
            % analyze data

            % loop over sessions
            for jj = 1:numSess
                % session index
                sessIndex = num2str(jj);
                % loop over components
                for kk = 1:numComp

                    % component index
                    fileIndex = icatb_returnFileIndex(kk);
                    if jj == 1
                        %%%Mean for different sessions and different subjects %%%
                        outputFiles(1).ses(1).name(kk,:) = [prefix, meanAllNaming, COMPONENT_NAMING, ...
                            SESSION_POSTFIX, '_all', '_', fileIndex, fileExtn];
                    end
                    % single subject single session component naming
                    outputFiles(2).ses(jj).name(kk, :) = ...
                        [prefix, subNaming, '01', COMPONENT_NAMING, SESSION_POSTFIX, sessIndex, ...
                        '_', fileIndex, fileExtn];

                end
                % end loop over components

            end
            % end loop over sessions

        elseif strcmpi(fileExtn, '.nii') || strcmpi(fileExtn, '.mat')
            % nifti data

            % loop over sessions
            for jj = 1:numSess
                % session index
                sessIndex = num2str(jj);

                if jj == 1
                    %%%Mean for different sessions and different subjects %%%
                    outputFiles(1).ses(1).name = [prefix, meanAllNaming, COMPONENT_NAMING, ...
                        SESSION_POSTFIX, '_all', '_', fileExtn];
                end
                % single subject single session component naming
                outputFiles(2).ses(jj).name = ...
                    [prefix, subNaming, '01', COMPONENT_NAMING, SESSION_POSTFIX, ...
                    sessIndex, '_', fileExtn];

            end
            % end loop over sessions

        else
            error('Unknown image extensions');
        end
    else
        % multiple subjects and multiple sessions

        % check image extensions
        if strcmpi(fileExtn, '.img')
            % analyze data

            % loop over subjects
            for ii = 1:numSub

                % subject index
                subIndex = icatb_returnFileIndex(ii);

                % loop over sessions
                for jj = 1:numSess
                    % session index
                    sessIndex = num2str(jj);
                    % loop over components
                    for kk = 1:numComp

                        % component index
                        fileIndex = icatb_returnFileIndex(kk);
                        if ii == 1
                            outputFiles(MEAN_INDEX).ses(jj).name(kk, :) = [prefix, meanNaming, COMPONENT_NAMING, ...
                                SESSION_POSTFIX, sessIndex, '_', fileIndex, fileExtn];
                            outputFiles(STD_INDEX).ses(jj).name(kk, :) = [prefix, stdNaming, COMPONENT_NAMING, ...
                                SESSION_POSTFIX, sessIndex, '_', fileIndex, fileExtn];
                            outputFiles(TMAP_INDEX).ses(jj).name(kk, :) = [prefix, tmapNaming, COMPONENT_NAMING, ...
                                SESSION_POSTFIX, sessIndex, '_', fileIndex, fileExtn];
                        end
                        if ii == 1 && jj == 1
                            outputFiles(MEAN_ALL_INDEX).ses(1).name(kk, :) = [prefix, meanAllNaming, COMPONENT_NAMING, ...
                                SESSION_POSTFIX, '_all', '_', fileIndex, fileExtn];
                        end
                        outputFiles(SUBJECT_ICA_INDEX + ii - 1).ses(jj).name(kk, :) = [prefix, subNaming,...
                            subIndex, COMPONENT_NAMING, SESSION_POSTFIX, sessIndex, '_', fileIndex, fileExtn];

                    end
                    % end loop over components

                end
                % end loop over sessions

            end
            % end loop over subjects

        elseif strcmpi(fileExtn, '.nii') || strcmpi(fileExtn, '.mat')
            % nifti data

            % loop over subjects
            for ii = 1:numSub

                % subject index
                subIndex = icatb_returnFileIndex(ii);

                % loop over sessions
                for jj = 1:numSess
                    % session index
                    sessIndex = num2str(jj);
                    if ii == 1
                        outputFiles(MEAN_INDEX).ses(jj).name = [prefix, meanNaming, COMPONENT_NAMING, ...
                            SESSION_POSTFIX, sessIndex, '_', fileExtn];
                        outputFiles(STD_INDEX).ses(jj).name = [prefix, stdNaming, COMPONENT_NAMING, ...
                            SESSION_POSTFIX, sessIndex, '_', fileExtn];
                        outputFiles(TMAP_INDEX).ses(jj).name = [prefix, tmapNaming, COMPONENT_NAMING, ...
                            SESSION_POSTFIX, sessIndex, '_', fileExtn];
                    end
                    if ii == 1 && jj == 1
                        outputFiles(MEAN_ALL_INDEX).ses(1).name = [prefix, meanAllNaming, COMPONENT_NAMING, ...
                            SESSION_POSTFIX, '_all', '_', fileExtn];
                    end
                    outputFiles(SUBJECT_ICA_INDEX + ii - 1).ses(jj).name = [prefix, subNaming,...
                        subIndex, COMPONENT_NAMING, SESSION_POSTFIX, sessIndex, '_', fileExtn];

                end
                % end loop over sessions

            end
            % end loop over subjects

        else
            error('Unknown image extensions');
        end

    end
    % end for checking data sets

else

    cNaming = [prefix, subNaming, COMPONENT_NAMING];

    if (strcmpi(fileExtn, '.nii'))
        outputFiles(1).ses(1).name = [cNaming, fileExtn];
    else
        files = cell(numComp, 1);
        for nComp = 1:numComp
            files{nComp} = [cNaming, icatb_returnFileIndex(nComp), fileExtn];
        end
        outputFiles(1).ses(1).name = char(files);
    end

end