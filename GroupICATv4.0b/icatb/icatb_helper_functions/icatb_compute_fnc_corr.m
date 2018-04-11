function fncVals = icatb_compute_fnc_corr(param_file, TRs, varargin)
%% Compute fnc correlations
%
% Inputs:
% 1. param_file - ICA parameter file (*ica*param*mat)
% 2. TRs - TR in seconds. If subjects have different TRs, enter TRs in a
% vector of length equal to the number of subjects.
% 3. varargin - Variable arguments passed in pairs:
%   a. filter_params - Filter cutoff in Hz. By default, low pass filter is used. If you specify two element vector,
%   bandpass filter is used.
%   b. lag - Lag in seconds. By default, fnc correlations are computed
%   using no lag.
%   c. comps - Component numbers in a vector.
%   d. shift_resolution - Shift resolution
%   e. despike - Options are 1 and 0.
%   f. subjects - Subject numbers in a vector.
%
% Outputs:
% fnc correlations - subjects by correlation values (fisher z-transformed).
% If lagged correlation is used, third dimension contains correlations (fnc_corrs(:,:,1)) and
% lag (fnc_corrs(:, :, 2)) respectively.
%

icatb_defaults;
global TIMECOURSE_POSTPROCESS;
global MANCOVA_DEFAULTS;

filter_params = 0.15;
try
    filter_params = TIMECOURSE_POSTPROCESS.fnc.cutoff_frequency;
catch
end

dMaxLag = 0;
try
    dMaxLag = MANCOVA_DEFAULTS.fnc.lag;
catch
end

nShiftResolution = 25;
try
    nShiftResolution = MANCOVA_DEFAULTS.fnc.shift_resolution;
catch
end

despikeTC = 0;
detrend_no = 1;
covariate_files = '';
scansToInclude = [];

for n = 1:2:length(varargin)
    if (strcmpi(varargin{n}, 'filter_params'))
        filter_params = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'lag'))
        dMaxLag = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'comps'))
        comps = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'shift_resolution'))
        nShiftResolution = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'despike'))
        despikeTC = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'subjects'))
        subjects = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'detrend_no'))
        detrend_no = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'covariates'))
        covariate_files = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'scansToInclude'))
        scansToInclude = varargin{n + 1};
    end
end


%% Get params
if (ischar(despikeTC) && strcmpi(despikeTC, 'yes'))
    despikeTC = 1;
end

if (ischar(param_file))
    
    load(param_file);
    pathstr = fileparts(param_file);
    sesInfo.outputDir = pathstr;
    sesInfo.userInput.pwd = sesInfo.outputDir;
    
else
    
    sesInfo = param_file;
    
end

% diffTimePoints = sesInfo.diffTimePoints;
% if (length(TRs) == 1)
%     TRs = TRs*ones(1, length(diffTimePoints));
% end

if (~exist('comps', 'var'))
    comps = 1:sesInfo.numComp;
end

if (~exist('subjects', 'var'))
    subjects = (1:sesInfo.numOfSub);
end

if (length(TRs) == 1)
    TRs = TRs*ones(1, length(subjects));
end

minTR = min(TRs);
dInterpBy = floor((minTR / dMaxLag) *nShiftResolution);
anComb = nchoosek(1:length(comps), 2);
nCombCnt = size(anComb,1);
rovdInnerLoop = -nShiftResolution:nShiftResolution;
adCorrTemp = zeros(sesInfo.numOfSub*sesInfo.numOfSess, nCombCnt, size(rovdInnerLoop, 2) );



nV = 0;
subCount = 0;
for nSub = subjects
    subCount = subCount + 1;
    for nSess = 1:sesInfo.numOfSess
        nV = nV + 1;
        currentTR = TRs(subCount);
        disp(['Loading subject #', num2str(nSub), ' session # ', num2str(nSess), ' ...']);
        tc = icatb_loadComp(sesInfo, comps, 'subjects', nSub, 'sessions', nSess, 'vars_to_load', 'tc', 'detrend_no', detrend_no, ...
            'covariates', covariate_files, 'scansToInclude', scansToInclude);
        tc = squeeze(tc);
        
        if (currentTR ~= minTR)
            interpFactor = currentTR/minTR;
            tc = icatb_interp_data(tc, interpFactor);
        end
        
        if (min(filter_params) > 0)
            tc = icatb_filt_data(tc, minTR, filter_params);
        end
        
        if (despikeTC)
            tc = icatb_despike_tc(tc, minTR);
        end
        
        %if (dMaxLag == 0)
        
        c = icatb_corr(tc);
        c(1:size(c, 1) + 1:end) = 0;
        c = icatb_mat2vec(c);
        
        
        if (nV == 1)
            fnc_corrs = zeros(length(subjects)*sesInfo.numOfSess, nCombCnt);
        end
        
        fnc_corrs(nV, :) = c(:)';
        
        % else
        
        if (dMaxLag > 0)
            % lagged correlation (from fnc toolbox)
            
            newTC = icatb_interp_data(tc, dInterpBy);
            
            for i = 1:nCombCnt
                rovdT1 = squeeze(newTC(:, anComb(i, 1)));
                rovdT2 = squeeze(newTC(:, anComb(i, 2)));
                rmrovdT1 = rovdT1 - mean(rovdT1);
                rmrovdT2 = rovdT2 - mean(rovdT2);
                [cc, ll] = xcorr(rmrovdT1, rmrovdT2, 'coeff');
                adCorrTemp(nV, i, 1:length(rovdInnerLoop)) = cc(find(ll==0)-nShiftResolution:find(ll==0)+nShiftResolution);
            end
            
            
        end
        
        clear newTC;
        
    end
    
end


fnc_corrs = reshape(fnc_corrs, sesInfo.numOfSess, length(subjects), nCombCnt);
fnc_corrs = permute(fnc_corrs, [2, 1, 3]);
fnc_corrs = icatb_r_to_z(fnc_corrs);

fncVals.no_lag.values = fnc_corrs;

if (dMaxLag > 0)
    % Lagged correlations
    
    [adCorrPos_v, adCorrPos_l] = max(adCorrTemp,[],3);
    [adCorrAbs_v, adCorrAbs_l ] = signedmax(adCorrTemp,3);
    adCorrAbs_v = icatb_r_to_z(adCorrAbs_v);
    adCorrPos_v = icatb_r_to_z(adCorrPos_v);
    
    adCorrPos(:,:,1) = adCorrPos_v;  adCorrPos(:,:,2) = adCorrPos_l;
    adCorrPos(:,:,3) = (adCorrPos_l - 1 - nShiftResolution)/dInterpBy*minTR;
    
    adCorrAbs(:,:,1) = adCorrAbs_v;  adCorrAbs(:,:,2) = adCorrAbs_l;
    adCorrAbs(:,:,3) = (adCorrAbs_l - 1 - nShiftResolution)/dInterpBy*minTR;
    
    adCorrAbs = reshape(adCorrAbs, sesInfo.numOfSess, length(subjects), nCombCnt, 3);
    adCorrAbs = permute(adCorrAbs, [2, 1, 3, 4]);
    adCorrPos = reshape(adCorrPos, sesInfo.numOfSess, length(subjects), nCombCnt, 3);
    adCorrPos = permute(adCorrPos, [2, 1, 3, 4]);
    
    fncVals.lag.dMaxLag = dMaxLag;
    fncVals.lag.absCorr = adCorrAbs;
    fncVals.lag.posCorr = adCorrPos;
    
end

disp('Done');


function [smax, sind ] = signedmax(data,dim)

if nargin < 2
    dim = 1;
end

[maxp, posind] = max(data, [], dim);
[maxn, negind] = max(-data, [], dim);

smax = maxp;
nind = find(maxn-maxp > 0);
smax(nind) = -maxn(nind);

sind = posind;
sind(nind) = negind(nind);
