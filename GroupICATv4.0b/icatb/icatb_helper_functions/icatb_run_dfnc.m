function dfncInfo = icatb_run_dfnc(param_file)
%% Run dfnc

icatb_defaults;
global DETRENDNUMBER;

if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('title', 'Select dFNC Parameter File', 'typeEntity', 'file', 'typeSelection', ...
        'single', 'filter', '*dfnc.mat');
    drawnow;
end

if (isempty(param_file))
    error('Please select dFNC parameter file');
end

if (ischar(param_file))
    load(param_file);
    outputDir = fileparts(param_file);
    if (isempty(outputDir))
        outputDir = pwd;
    end
    if (~exist('dfncInfo', 'var'))
        error(['Selected file ', param_file, ' is not a valid dFNC parameter file']);
    end
else
    dfncInfo = param_file;
    outputDir = dfncInfo.userInput.outputDir;
    clear param_file;
end

cd (outputDir);

%% Initialise
outputDir = pwd;
detrend_no = DETRENDNUMBER;
doDespike = 'no';
tc_filter = 0;
method = 'none';
wsize = 30;
window_alpha = 3;
window_type = 'gaussian';
num_repetitions = 10;
covariates = 'none';

try
    detrend_no = dfncInfo.userInput.feature_params.final.tc_detrend;
    doDespike = dfncInfo.userInput.feature_params.final.tc_despike;
    tc_filter = dfncInfo.userInput.feature_params.final.tc_filter;
    wsize = dfncInfo.userInput.feature_params.final.wsize;
    window_alpha = dfncInfo.userInput.feature_params.final.window_alpha;
    num_repetitions = dfncInfo.userInput.feature_params.final.num_repetitions;
    method = dfncInfo.userInput.feature_params.final.method;
    covariates = dfncInfo.userInput.feature_params.final.tc_covariates;
    window_type = dfncInfo.userInput.feature_params.final.window_type;
catch
end

%% Tukey window option is disabled for now.
if (strcmpi(window_type, 'tukey'))
    %if (window_alpha <= 1)
    window_alpha = 3;
    %end
    disp(['Tukey window option is currently disabled. Using gaussian window with alpha = ', ...
        num2str(window_alpha), ' ...']);
end


TR = dfncInfo.userInput.TR;

load(dfncInfo.userInput.ica_param_file);

comps = [dfncInfo.userInput.comp.value];

if (length(comps) ~= length(unique(comps)))
    error('One or more components are replicated across different component networks');
end

if ((length(TR) > 1) && (length(TR) ~= sesInfo.numOfSub))
    error(['You have specified multiple TRs. TRs should match the length of number of subjects (', num2str(sesInfo.numOfSub), ')']);
end

dfncInfo.prefix = dfncInfo.userInput.prefix;
dfncInfo.comps = comps;
dfncInfo.outputDir = outputDir;
dfncInfo.TR = TR;
dfncInfo.detrend_no = detrend_no;
dfncInfo.doDespike = doDespike;
dfncInfo.tc_filter = tc_filter;
dfncInfo.method = method;
dfncInfo.wsize = wsize;
dfncInfo.window_alpha = window_alpha;
dfncInfo.window_type = window_type;
dfncInfo.num_repetitions = num_repetitions;

outDir = fileparts(dfncInfo.userInput.ica_param_file);
sesInfo.outputDir = outDir;


tic;
logFile = fullfile(outputDir, [dfncInfo.prefix, '_results.log']);
diary(logFile);

disp('----------------------------------------------------');
disp('-------------- STARTING DYNAMIC FNC --------------');
disp('----------------------------------------------------');

disp('Loading timecourses ....');

covariate_files = '';
scansToInclude = [];
if (~strcmpi(covariates, 'none'))
    try
        covariate_files = dfncInfo.userInput.feature_params.final.tc_covariates_userdata.filesList;
        scansToInclude = dfncInfo.userInput.feature_params.final.tc_covariates_userdata.file_numbers;
        disp('Variance associated with the covariates will be removed from the timecourses ... ');
    catch
    end
end

fprintf('\n');

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

%% Make sliding window
numComp = length(comps);
c = compute_sliding_window(nT, window_alpha, wsize);
A = repmat(c, 1, numComp);

Nwin = nT - wsize;
%initial_lambdas = (0.1:.03:1);
initial_lambdas = (0.1:.03:.40);
best_lambda = zeros(1, sesInfo.numOfSub);

outFiles = cell(1, sesInfo.numOfSub);

tc = reshape(tc, sesInfo.numOfSub, sesInfo.numOfSess, nT, numComp);

dataSetCount = 0;

%% Loop over subjects
for nSub = 1:sesInfo.numOfSub
    
    %% Loop over sessions
    for nSess = 1:sesInfo.numOfSess
        
        dataSetCount = dataSetCount + 1;
        
        results_file = [dfncInfo.userInput.prefix, '_sub_', icatb_returnFileIndex(nSub), '_sess_', icatb_returnFileIndex(nSess), '_results.mat'];
        
        FNCdyn = zeros(Nwin, numComp*(numComp - 1)/2);
        Lambdas = zeros(num_repetitions, length(initial_lambdas));
        disp(['Computing dynamic FNC on subject ', num2str(nSub), ' session ', num2str(nSess)]);
        
        if (strcmpi(doDespike, 'yes') && (tc_filter > 0))
            disp('Despiking and filtering timecourses ...');
        elseif (strcmpi(doDespike, 'yes') && (tc_filter == 0))
            disp('Despiking timecourses ...');
        elseif (strcmpi(doDespike, 'no') && (tc_filter > 0))
            disp('Filtering timecourses ...');
        end
        
        %% Preprocess timecourses
        for nComp = 1:numComp
            current_tc = squeeze(tc(nSub, nSess, :, nComp));
            
            % Despiking timecourses
            if (strcmpi(doDespike, 'yes'))
                current_tc = icatb_despike_tc(current_tc, TR);
            end
            
            % Filter timecourses
            if (tc_filter > 0)
                current_tc = icatb_filt_data(current_tc, TR, tc_filter);
            end
            tc(nSub, nSess, :, nComp) = current_tc;
        end
        
        %% Apply circular shift to timecourses
        tcwin = zeros(Nwin, nT, numComp);
        for ii = 1:Nwin
            Ashift = circshift(A, round(-nT/2) + round(wsize/2) + ii);
            tcwin(ii, :, :) = squeeze(tc(nSub, nSess, :, :)).*Ashift;
        end
        
        if strcmpi(method, 'L1')
            
            useMEX = 0;
            
            try
                GraphicalLassoPath([1, 0; 0, 1], 0.1);
                useMEX = 1;
            catch
            end
            
            disp('Using L1 regularisation ...');
            
            %% L1 regularisation
            Pdyn = zeros(Nwin, numComp*(numComp - 1)/2);
            
            fprintf('\t rep ')
            %% Loop over no of repetitions
            for r = 1:num_repetitions
                fprintf('%d, ', r)
                [trainTC, testTC] = split_timewindows(tcwin, 1);
                trainTC = icatb_zscore(trainTC);
                testTC = icatb_zscore(testTC);
                [wList, thetaList] = computeGlasso(trainTC, initial_lambdas, useMEX);
                obs_cov = icatb_cov(testTC);
                L = cov_likelihood(obs_cov, thetaList);
                Lambdas(r, :) = L;
            end
            
            fprintf('\n')
            [mv, minIND] =min(Lambdas, [], 2);
            blambda = mean(initial_lambdas(minIND));
            fprintf('\tBest Lambda: %0.3f\n', blambda)
            best_lambda(dataSetCount) = blambda;
            
            % now actually compute the covariance matrix
            fprintf('\tWorking on estimating covariance matrix for each time window...\n')
            for ii = 1:Nwin
                %fprintf('\tWorking on window %d of %d\n', ii, Nwin)
                [wList, thetaList] = computeGlasso(icatb_zscore(squeeze(tcwin(ii, :, :))), blambda, useMEX);
                a = icatb_corrcov(wList);
                a = a - eye(numComp);
                FNCdyn(ii, :) = mat2vec(a);
                InvC = -thetaList;
                r = (InvC ./ repmat(sqrt(abs(diag(InvC))), 1, numComp)) ./ repmat(sqrt(abs(diag(InvC)))', numComp, 1);
                r = r + eye(numComp);
                Pdyn(ii, :) = mat2vec(r);
            end
            
            FNCdyn = atanh(FNCdyn);
            disp(['.... saving file ', results_file]);
            icatb_save(fullfile(outputDir, results_file), 'Pdyn', 'Lambdas', 'FNCdyn');
            
        elseif strcmpi(method, 'none')
            %% No L1
            for ii = 1:Nwin
                a = icatb_corr(squeeze(tcwin(ii, :, :)));
                FNCdyn(ii, :) = mat2vec(a);
            end
            FNCdyn = atanh(FNCdyn);
            
            disp(['.... saving file ', results_file]);
            icatb_save(fullfile(outputDir, results_file), 'FNCdyn');
        end
        
        outFiles{dataSetCount} = results_file;
        disp('Done');
        fprintf('\n');
        
    end
    %% End of loop over sessions
end
%% End of loop over subjects

dfncInfo.best_lambda = best_lambda;
dfncInfo.outputFiles = outFiles;
fileN = fullfile(outputDir, [dfncInfo.prefix, '.mat']);
disp(['Saving parameter file ', fileN, ' ...']);
icatb_save(fileN, 'dfncInfo');

totalTime = toc;

fprintf('\n');

disp('Analysis Complete');

disp(['Total time taken to complete the analysis is ', num2str(totalTime/60), ' minutes']);

diary('off');

fprintf('\n');

function c = compute_sliding_window(nT, win_alpha, wsize)
%% Compute sliding window
%

nT1 = nT;
if mod(nT, 2) ~= 0
    nT = nT + 1;
end

m = nT/2;
w = round(wsize/2);
%if (strcmpi(win_type, 'tukey'))
%    gw = icatb_tukeywin(nT, win_alpha);
%else
gw = gaussianwindow(nT, m, win_alpha);
%end
b = zeros(nT, 1);  b((m -w + 1):(m+w)) = 1;
c = conv(gw, b); c = c/max(c); c = c(m+1:end-m+1);
c = c(1:nT1);


function [vec, IND] = mat2vec(mat)
% vec = mat2vec(mat)
% returns the lower triangle of mat
% mat should be square

[n,m] = size(mat);

if n ~=m
    error('mat must be square!')
end


temp = ones(n);
%% find the indices of the lower triangle of the matrix
IND = find((temp-triu(temp))>0);

vec = mat(IND);


% function w = gaussianwindow(N,x0,sigma)
%
% x = 0:N-1;
% w = exp(- ((x-x0).^2)/ (2 * sigma * sigma))';


function L = cov_likelihood(obs_cov, theta)
% L = cov_likelihood(obs_cov, sigma)
% INPUT:
% obs_cov is the observed covariance matrix
% theta is the model precision matrix (inverse covariance matrix)
% theta can be [N x N x p], where p lambdas were used
% OUTPUT:
% L is the negative log-likelihood of observing the data under the model
% which we would like to minimize

nmodels = size(theta,3);

L = zeros(1,nmodels);
for ii = 1:nmodels
    % log likelihood
    theta_ii = squeeze(theta(:,:,ii));
    L(ii) = -log(det(theta_ii)) + trace(theta_ii*obs_cov);
end

function [trainTC, testTC] = split_timewindows(TCwin, ntrain)
%[Nwin, nT, nC] = size(TCwin);


[Nwin, nT, nC] = size(TCwin);

r = randperm(Nwin);
trainTC = TCwin(r(1:ntrain),:,:);
testTC = TCwin(r(ntrain+1:end),:,:);

trainTC = reshape(trainTC, ntrain*nT, nC);
testTC = reshape(testTC, (Nwin-ntrain)*nT, nC);


function w = gaussianwindow(N,x0,sigma)

x = 0:N-1;
w = exp(- ((x-x0).^2)/ (2 * sigma * sigma))';


function [wList, thetaList] = computeGlasso(tc, initial_lambdas, useMEX)
%% Compute graphical lasso


if (useMEX == 1)
    [wList, thetaList] = GraphicalLassoPath(tc, initial_lambdas);
else
    tol = 1e-4;
    maxIter = 1e4;
    S = icatb_cov(tc);
    thetaList = zeros(size(S, 1), size(S, 2), length(initial_lambdas));
    wList = thetaList;
    
    for nL = 1:size(wList, 3)
        [thetaList(:, :, nL), wList(:, :, nL)] = icatb_graphicalLasso(S, initial_lambdas(nL), maxIter, tol);
    end
    
end


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