function [tc, TR] = icatb_loadAndInterpTC(sesInfo, comps, TR, varargin)
%% Load and interpolate TCs for dFNC computation
%

icatb_defaults;
global DETRENDNUMBER;

detrend_no = DETRENDNUMBER;
covariate_files = '';
scansToInclude = [];

for i = 1:2:length(varargin)
    if (strcmpi(varargin{i}, 'detrend_no'))
        detrend_no = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'covariates'))
        covariate_files = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'scansToInclude'))
        scansToInclude = varargin{i + 1};
    end
end


%% Load timecourses
if ((length(TR) == 1) || all(TR == min(TR)))
    nT = min(sesInfo.diffTimePoints);
    if (all(sesInfo.diffTimePoints  == nT))
        % Load timecourses
        tc = icatb_loadComp(sesInfo, comps, 'vars_to_load', 'tc', 'detrend_no', detrend_no, 'truncate_tp', 1, 'covariates', ...
            covariate_files, 'scansToInclude', scansToInclude);
        
    else
        
        nT = min(sesInfo.diffTimePoints);
        tcs = cell(sesInfo.numOfSub, sesInfo.numOfSess);
        %% Loop over subjects
        for nSub = 1:sesInfo.numOfSub
            %% Loop over sessions
            for nSess = 1:sesInfo.numOfSess
                tmp = icatb_loadComp(sesInfo, comps, 'subjects', nSub, 'sessions', nSess, 'vars_to_load', 'tc', 'detrend_no', detrend_no, 'covariates', ...
                    covariate_files, 'scansToInclude', scansToInclude);
                tmp = squeeze(tmp);
                tcs{nSub, nSess} = tmp;
            end
        end
        tc = getTruncatedTcs(tcs, nT);
        clear tcs;
    end
    
else
    % Interpolate timecourses using minimum TR and maximum TR. Truncate the
    % timepoints to match the same between subjects and sessions
    disp('Interpolating timecourses using the minimum TR and the dataset TR ...');
    tcs = cell(sesInfo.numOfSub, sesInfo.numOfSess);
    
    countTmp = 0;
    %% Loop over subjects
    for nSub = 1:sesInfo.numOfSub
        %% Loop over sessions
        for nSess = 1:sesInfo.numOfSess
            countTmp = countTmp + 1;
            tmp = icatb_loadComp(sesInfo, comps, 'subjects', nSub, 'sessions', nSess, 'vars_to_load', 'tc', 'detrend_no', detrend_no, 'covariates', ...
                covariate_files, 'scansToInclude', scansToInclude);
            tmp = squeeze(tmp);
            interpFactor = TR(nSub)/min(TR);
            [num, denom] = rat(interpFactor);
            tmp = resample(tmp, num, denom);
            tmpLen = size(tmp, 1);
            if (countTmp == 1)
                nT = tmpLen;
            else
                nT = min([nT, tmpLen]);
            end
            tcs{nSub, nSess} = tmp;
        end
        %% End of Loop over sessions
    end
    %% End of Loop over subjects
    
    
    %% Truncate timeseries to match subjects and sessions
    tc = getTruncatedTcs(tcs, nT);
    %     tc = zeros(sesInfo.numOfSub, sesInfo.numOfSess, nT, length(comps));
    %     countTmp = 0;
    %     %% Loop over subjects
    %     for nSub = 1:sesInfo.numOfSub
    %         %% Loop over sessions
    %         for nSess = 1:sesInfo.numOfSess
    %             countTmp = countTmp + 1;
    %             tmp = tcs{nSub, nSess};
    %             tmp = tmp(1:nT, :);
    %             tc(nSub, nSess, :, :) = tmp;
    %         end
    %     end
    clear tcs;
end

TR = min(TR);


function tc = getTruncatedTcs(tcs, nT)


numOfSub = size(tcs, 1);
numOfSess = size(tcs, 2);
numComp = size(tcs{1,1}, 2);
tc = zeros(numOfSub, numOfSess, nT, numComp);

%% Loop over subjects
for nSub = 1:numOfSub
    %% Loop over sessions
    for nSess = 1:numOfSess
        tmp = tcs{nSub, nSess};
        tmp = squeeze(tmp);
        tp = size(tmp, 1);
        if (tp ~= nT)
            interpFactor = nT/tp;
            [num, denom] = rat(interpFactor);
            tmp = resample(tmp, num, denom);
        end
        tc(nSub, nSess, :, :) = tmp;
    end
end