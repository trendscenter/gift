function dataStruct = icatb_check_pca_opts(dataStruct)
%% Remove covariance opts field
%
if (isfield(dataStruct, 'pca_opts'))
    pca_opts = dataStruct.pca_opts;
else
    if (isfield(dataStruct, 'covariance_opts'))
        pca_opts = dataStruct.covariance_opts;
        dataStruct = rmfield(dataStruct, 'covariance_opts');
    end
end

if (~exist('pca_opts', 'var'))
    pca_opts = icatb_pca_options(dataStruct.pcaType);
else
    pca_opts = icatb_pca_options(dataStruct.pcaType, pca_opts, 'off');
end


if (isnumeric(dataStruct.pcaType))
    tmp = icatb_pca_options;
    tmp = cellstr(tmp);
    dataStruct.pcaType = lower(tmp{dataStruct.pcaType});
end

pca_opts.precision = icatb_checkPrecision(pca_opts.precision);

dataStruct.pca_opts = pca_opts;