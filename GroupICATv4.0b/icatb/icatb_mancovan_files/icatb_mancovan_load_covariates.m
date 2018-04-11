function val = icatb_mancovan_load_covariates(txtFile, covType, numOfSub)
%% Load covariates
%

txtFile = deblank(txtFile);

if (strcmpi(covType, 'continuous'))
    val = load(txtFile, '-ascii');
    if (prod(size(val)) ~= length(val))
        error('The selected file contains a matrix. Only vector is allowed');
    end
    val = val(:);
    val = num2str(val);
else
    fid = fopen(txtFile, 'r');
    if (fid == -1)
        error(['File ', txtFile, ' cannot be opened for reading']);
    end
    try
        dd = textscan(fid, '%s', 'delimiter', '\t\n,', 'multipleDelimsAsOne', 1);
        val = dd{1};
    catch
        val = [];
    end
    fclose(fid);
end

val = strtrim(cellstr(val));

val = val(:)';

if (exist('numOfSub', 'var'))
    if (length(val) ~= numOfSub)
        error('Covariate vector must match the no. of subjects in the analysis');
    end
end