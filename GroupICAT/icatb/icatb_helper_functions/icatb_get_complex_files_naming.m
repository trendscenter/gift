function [P] = icatb_get_complex_files_naming(P, dataType, complexInfo, optional)
% name complex images
%
% Input:
% 1. P - files (of string type)
% 2. dataType - real or complex
% 3. complexInfo - contains two fields like complexType and complexNaming
% 4. optional - files are read or written
%
% Ouput:
% P - string for real images and a structure with fields first and second
% for complex images.

icatb_defaults;
% FUNTIONAL FILE DEFAULTS
global FUNCTIONAL_DATA_FILTER;

if ~exist('dataType', 'var')
    dataType = 'real';
end

if ~exist('optional', 'var')
    optional = 'read';
end


% if the data type is complex then form real and imaginary or magnitude and
% phase namings
if strcmpi(dataType, 'complex')

    % complex file naming
    complex_file_naming = complexInfo.complex_file_naming;

    % replicate the structure
    P1.first = [];
    P1.second = [];
    P1 = repmat(P1, size(P, 1), 1);
    % loop over files
    for ii = 1:size(P, 1)
        % get the file namings
        [pathstr, name, extn] = fileparts(deblank(P(ii, :)));
        % check for underscore positions
        underScorePositions = icatb_findstr(name, '_');
        % put a check for underscore
        if ~isempty(underScorePositions)
            % if the files are of read type
            if strcmpi(optional, 'read')

                % compare the prefix with the second file naming in the defaults
                if strcmp(name(1:underScorePositions(1)), complex_file_naming{2})
                    % form prefixes for only first part of complex images
                    P1(ii).second = fullfile(pathstr, [name, extn]);
                    P1(ii).first = fullfile(pathstr, [complex_file_naming{1}, name(underScorePositions(1)+1:end), ...
                        extn]);

                elseif strcmp(name(underScorePositions(end):end), complex_file_naming{2})

                    % form suffixes for the first part
                    P1(ii).second = fullfile(pathstr, [name, extn]);
                    P1(ii).first = fullfile(pathstr, [name(1:underScorePositions(end)-1), ...
                        complex_file_naming{1}, extn]);

                else
                    error('please check the complex file naming in icatb_defaults.m');
                end
                % end for checking cases for read type
            end
            % end for checking read type

            % if the files are of read type
            if strcmpi(optional, 'write')
                
                if isempty(extn)
                    extn = FUNCTIONAL_DATA_FILTER;
                    if strcmp(extn(1), '*')
                        extn = extn(2:end);
                    end

                end

                % form prefixes for both the first and second
                P1(ii).second = fullfile(pathstr, [complex_file_naming{2}, name, extn]);
                P1(ii).first = fullfile(pathstr, [complex_file_naming{1}, name, extn]);
                % % compare the suffix with the second file naming in the
                % defaults
            end
            % end for checking the write type

        else
            error('file names should contain underscore to distinguish complex images');
        end
        % end for checking underscore positions

    end
    % end loop over files
    P = P1;

    clear P1;
end
% end for checking data type