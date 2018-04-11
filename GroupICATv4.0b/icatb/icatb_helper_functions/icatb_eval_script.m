function inputData = icatb_eval_script(file_name)
%% Evaluate script and store the variables in a data structure
%
% Inputs:
% file_name - File name
%
% Outputs:
% inputData - Input data structure
%

if (isdeployed)
    % execute scripts differently in deployed mode
    fid = fopen(file_name, 'r');
    
    if (fid == -1)
        error(['File ', file_name, ' cannot be opened']);
    end
    
    try
        tmp_strs = fread(fid, '*char');
        fclose(fid);
    catch
        fclose(fid);
    end
    
    eval(tmp_strs');
    
    clear fid tmp_strs;
    
else
    
    oldDir = pwd;
    
    %% Do file parts
    [pathstr, fName, extn] = fileparts(file_name);
    if isempty(pathstr)
        pathstr = pwd;
    end
    
    cd(pathstr);
    
    %% Evaluate file
    eval(fName);
    
    cd(oldDir);
    
end

if (nargout == 1)
    
    vars = whos;
    
    if (length(vars) == 1)
        error(['No variables found in file ', file_name]);
    end
    
    % Generate inputData
    inputData = struct;
    for n = 1:length(vars)
        inputData = setfield(inputData, vars(n).name, eval(vars(n).name));
    end
    
end