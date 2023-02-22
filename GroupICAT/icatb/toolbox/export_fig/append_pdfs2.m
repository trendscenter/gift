function append_pdfs2(varargin)
% APPEND_PDFS2 Appends/concatenates multiple PDF files
% NOTE! 041622 version of old append_pdfs.m, improving Windows compatability!!
%
% Example:
%   append_pdfs(outputFilename, inputFilename1, inputFilename2, ...)
%   append_pdfs(outputFilename, inputFilenames_list{:})
%   append_pdfs(outputFilename, inputFilenames_cell_or_string_array)
%   append_pdfs output.pdf input1.pdf input2.pdf
%
% This function appends multiple PDF files to an existing PDF file, or
% concatenates them into a PDF file if the output file doesn't yet exist.
%
% This function concatenates multiple PDF files to an existing PDF file
%
% This function only requires MATLAB
%
% IN:
%    output - string of output file name (including the extension, .pdf).
%    input_list - cell array list of input file name strings. All input
%                 files are appended in order.
% OUT:
%    On disk the file out.pdf is created
%
% This function originates from MATLAB File Exchange 
% Benjamin Groﬂmann (2022). Merge PDF-Documents 
% (https://www.mathworks.com/matlabcentral/fileexchange/89127-merge-pdf-documents), 
% MATLAB Central File Exchange. Retrieved April 15, 2022.

if nargin < 2,  return;  end  % sanity check

% Convert strings => chars; strtrim extra spaces
    varargin = cellfun(@str2char,varargin,'un',false);

    % Convert cell array into individual strings (issue #299)
    if nargin==2 && iscell(varargin{2})
        varargin = {varargin{1} varargin{2}{:}}; %#ok<CCAT>
    end

    % Are we appending or creating a new file
    append = exist(varargin{1}, 'file') == 2;
    output = [tempname '.pdf'];
    try
        % Ensure that the temp dir is writable (Javier Paredes 26/2/15)
        fid = fopen(output,'w');
        fwrite(fid,1);
        fclose(fid);
        delete(output);
        isTempDirOk = true;
    catch
        % Temp dir is not writable, so use the output folder
        [dummy,fname,fext] = fileparts(output); %#ok<ASGLU>
        fpath = fileparts(varargin{1});
        output = fullfile(fpath,[fname fext]);
        isTempDirOk = false;
    end
    if ~append
        output = varargin{1};
        varargin = varargin(2:end);
    end

    % part from files exchenge to merge pdfs
    memSet = org.apache.pdfbox.io.MemoryUsageSetting.setupMainMemoryOnly();
    merger = org.apache.pdfbox.multipdf.PDFMergerUtility;
    cellfun(@(f) merger.addSource(f), varargin)
    merger.setDestinationFileName(output)
    merger.mergeDocuments(memSet)

    % Rename the file if needed
    if append
        movefile(output, varargin{1}, 'f');
    end

end

% Convert a possible string => char
function value = str2char(value)
    try
        value = controllib.internal.util.hString2Char(value);
    catch
        if isa(value,'string')
            value = char(value);
        end
    end
    value = strtrim(value);
end