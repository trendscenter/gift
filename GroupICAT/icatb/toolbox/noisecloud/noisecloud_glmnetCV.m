function RESULT = noisecloud_glmnetCV(X,y,K,Lambda,alphav,type)
% noisecloud_glmnetCV for noisecloud package
% performs cross validated elastic net logistic regression

% X : data;
% y: class labels;
% K = number of cross validation folds (usually K=10)
% Lambda = [0.0625 0.125 0.25 0.5 1 2 4 8 16];
% alphav = 1:-0.01:0.1;
% type = 'binomial' or 'multinomial'

% Alpha which varies from 0 to 1 controls the sparsity. Alpha = 1 shall
% give you the sparsest solution i.e. most betas that are close to zero
% are made zero. Alpha = 0 gives you a simple ridge regression solution
% with no sparsification i.e betas of correlated variable are made
% equal. Unlike alpha, lambda doesn't have an intuitive correlate and is
% mostly a mathematical trick to control sparse and ridge regression at
% the same time.

% To know more, you can read Section 2 of
% http://www.jstatsoft.org/v33/i01/paper

options = glmnetSet;     % Initialize an empty glmnet options data structure
options.alpha = alphav;  % Fill it with our alpha values input

N = length(y);           % Number of Observations (in this case, networks)
[Indx,leftOver] = splitData(N,K);  % Split the N data observations into K sets

I = 1:K;

% If the user has not specified a set of values for Lambda, we will use
% glmnet to estimate a range.  If the user has specified Lambda, we set
% this in our options alpha.  The first takes much more time, so it might
% make sense to not specify it to get an ideal of the best values, and then
% to set those values in advance.

% If Lambda is specified
if ~isempty(Lambda)
    Lambda = fliplr(Lambda);
    opts.lambda = Lambda;
    opts.nlambda = length(Lambda);
    opts.lambda_min = min(Lambda);
    options = glmnetSet(opts);
    
    % If no lambda is specified (takes long time)
else
    fit = glmnet(X,y,type,options);
    Lambda = fit.lambda;
    opts.lambda = Lambda;
    opts.nlambda = length(Lambda);
    opts.lambda_min = min(Lambda);
    options = glmnetSet(opts);
end

% We we keep track of cva across all lambas and alphas to find the best
% combination, in a grid search
cum_cva = zeros(length(alphav),length(Lambda));

% For each fold, break apart the test and training set
fprintf('%s\n','Performing cross validated grid search for alpha...');
for fold = 1:K
    fprintf('%s\n',[ '     Running fold ' num2str(fold) ' of ' num2str(K) ]);
    ix = find(I ~= fold);  % Group (K) indices of "everyone else," training group
    trIx = Indx(ix,:); % Train Indices, produces an N-1 by M size matrix, each row is
    % an entire M subset, when we do trIX(:) in the next
    % line, we make into one column vector
    trIx = trIx(:);
    trIx = [trIx;leftOver'];  % We then add the leftover to the training group
    tstIx = Indx(fold,:);     % The test group is the original indexed data at the fold)
    tstIx = tstIx(:);         % turned into a column vector
    
    % Now, apply the indices to select the actual data subsets for test and train!
    Xtrain = X(trIx,:); % Training X
    ytrain = y(trIx);   % Training y
    
    % Now we do a grid search for alpha.
    for a = 1:length(alphav)
        options.alpha = alphav(a);                % Select alpha value
        
        
        fit = glmnet(Xtrain,ytrain,type,options); % Fit model using alpha, train and test specified above
        Xtest = X(tstIx,:);                       % Get the test data and labels
        ytest = y(tstIx);
        Yest = glmnetPredict(fit, 'class', Xtest); % predict the response, result Yest is a M x nLambas matrix
        Ytest = repmat(ytest,1,length(fit.lambda));    % each column is all predicted subject responses for particular Lambda
        Indxc = (Ytest == Yest);                     % Indxc will = 1 where test = predicted, 0 where not.  Each column is a particular lambda
        cva = sum(Indxc,1)./length(ytest);         % cross validation accuracy for particular alpha, a, across all lambdas
        cum_cva(a,1:length(cva)) = cum_cva(a,1:length(cva)) + cva; % Accumulate squared error for each fold
    end
end

cum_cva = cum_cva./K;
% Find Optimal Parameters
no_lambda = length(fit.lambda);
best_cva = 0;
for a = 1:length(alphav)
    for l = 1:no_lambda
        CVA = cum_cva(a,l);
        if CVA > best_cva
            opt_alpha = alphav(a);
            opt_lambda = Lambda(l);
            best_cva = CVA;
        end
    end
end

% Train with Optimal Parameters
fprintf('\n%s\n','Training cross validated model with optimal parameters...');
fprintf('%s\n',[ '     Alpha: ' num2str(opt_alpha) ]);
fprintf('%s\n',[ '     Best CVA: ' num2str(best_cva) ]);
fprintf('%s\n\n',[ '     Lambda: ' num2str(opt_lambda) ]);

options.alpha = opt_alpha;
fit_opt = glmnet(X,y,type,options);
wopt = glmnetCoef(fit_opt,opt_lambda);

% Save optimal parameters to result structure
RESULT.fit_opt = fit_opt;
RESULT.wopt = wopt;
RESULT.opt_alpha = opt_alpha;
RESULT.best_cva = best_cva;
RESULT.cum_cva = cum_cva;
RESULT.opt_lambda = opt_lambda;

% ROC CURVE WITH OPTIMAL LAMBDA AND ALPHA
% Now that we have cross validated labels and optimal values, plot the ROC curve
% varying across different "thresholds"

% Create data structure to hold predictions for each fold from other folds
Ypredictions = zeros(length(y),1);

% Set our options to only include the optimal lambda
options.nlambda = 1;
options.lambda = opt_lambda;
options.lambda_min = opt_lambda;

% For each fold, break apart the test and training set
for fold = 1:K
    ix = find(I ~= fold);  % Group (K) indices of "everyone else," training group
    trIx = Indx(ix,:);     % Train Indices, produces an N-1 by M size matrix, each row is
    % an entire M subset, when we do trIX(:) in the next
    % line, we make into one column vector
    trIx = trIx(:);
    trIx = [trIx;leftOver'];  % We then add the leftover to the training group
    tstIx = Indx(fold,:);     % The test group is the original indexed data at the fold)
    tstIx = tstIx(:);         % turned into a column vector
    
    % Now, apply the indices to select the actual data subsets for test and train!
    Xtrain = X(trIx,:); % Training X
    ytrain = y(trIx);   % Training y
    
    fit = glmnet(Xtrain,ytrain,type,options); % Fit model using alpha, train and test specified above
    
    Xtest = X(tstIx,:);                       % Get the test data and labels
    ytest = y(tstIx);
    Yres = glmnetPredict(fit, 'response', Xtest); % predict the response, result Yest is a M x nLambas matrix
    % Yres holds the predicted response, from 0 to 1. Closer to 0
    % corresponds to closer to label 1, or BAD.  Closer to 1 corresponds to
    % 2, or GOOD.  We will save this to our Ypredictions vector to create a
    % ROC curve after we finish with each fold.
    Ypredictions(tstIx) = Yres;
    
end

fprintf('%s\n','Plotting ROC curve');

% This figure will plot our ROC curve
f = figure;
subplot(1,2,1);
axis([0 1 0 1]); hold on
xlabel('1-specificity');
ylabel('sensitivity');

% Convert our correct labels to 1 for bad, 0 for good (instead of 2)
% Since we break into 10 and have some leftover, we will not use them to
% make the roc curve, so we set the number of evallabels = the number of
% predicted that we have
evallabels = y(1:length(Ypredictions));
evallabels(evallabels==2) = 0;

% To hold sensitivity and specificity values to calculate AROC
aroc_sens = [];
aroc_spec = [];

% Now set different thresholds, from 0(BAD) to 1(GOOD), to come up with "final"
% labels.  This is what we will vary for our ROC curve.
for i=0:.005:1
    % Set thresholds so 1=bad, 2=good (our original labels)
    Ypredicty = (Ypredictions < i);  % Now bad is 1 and 0 is good
    % Evaluate(ACTUAL,PREDICTED) returns EVAL = [accuracy sensitivity specificity precision recall f_measure gmean tp tn fp fn]
    EVAL = Evaluate(evallabels,Ypredicty);
    % Eval at 3 is specificity, Eval at 2 is sensitivity
    % Save to a matrix to calculate area under ROC
    aroc_sens = [ aroc_sens EVAL(2) ];
    aroc_spec = [ aroc_spec 1- EVAL(3) ];
end

% Now plot Aroc Sens, Aroc Spec
plot(aroc_spec,aroc_sens,'LineWidth',2);

% Calculate the area under the cross validated curve:
RESULT.auc = trapz(aroc_spec,aroc_sens);

% If you want to check the ROC curve, use MATLAB's function (only for one
% this is ONLY for one set of predictions, is NOT the cross validated ROC
% curve) - so the curves will be different
% [mX,mY,mT,mAUC,mOPTROCPT,mSUBY,mSUBYNAMES] = perfcurve(evallabels,Ypredictions,0);

% Add the "by chance" line
line(0:1,0:1,'color','r')
title([ 'ROC Curve for alpha= ' num2str(RESULT.opt_alpha) ', lambda= ' num2str(RESULT.opt_lambda) ]);

% Create a confusion matrix to go with best classifier.  We will set
% threshold to be .5, which is standard
Ypredicty = (Ypredictions < .5);  % Now bad is 1 and 0 is good
% Evaluate(ACTUAL,PREDICTED) returns EVAL = [accuracy sensitivity specificity precision recall f_measure gmean tp tn fp fn]
EVAL = Evaluate(evallabels,Ypredicty);
RESULT.ROC.acc = EVAL(1);
RESULT.ROC.sens = EVAL(2);
RESULT.ROC.spec = EVAL(3);
RESULT.ROC.tp = EVAL(8);
RESULT.ROC.tn = EVAL(9);
RESULT.ROC.fp = EVAL(10);
RESULT.ROC.fn = EVAL(11);

% Confusion Matrix
%f = figure('Position',[100 100 300 150]);
% data = [TP FP; FN TN]
subplot(1,2,2);
data = [EVAL(8) EVAL(10); EVAL(11) EVAL(9) ];
cnames = {'YES','NO'};
rnames = {'YES','NO'};
set(gca,'Units','normalized');
t = uitable('Data',data,'ColumnName',cnames,'RowName',rnames,'Tag','Table 4','TooltipString','Table 4','Parent',f,'Units','Normalized','Position',[.5 .5 .5 .3]);
title('Confusion Matrix'); axis off

% This function randomly splits N data into K groups of size M, and returns
% an ordered column vector of the indices, plus the "leftover" indices
function [indx,leftOver] = splitData(N,K)

M = floor(N/K);            % Number of Samples in each group
try
    Nr = randsample(N,N)';     % Create a row vector of length N with random # from
catch
    Nr = randperm(N);
    Nr = Nr(:)';
end
indx = [];                 % 1 to N assigned to each spot... this is our random sample
cnt = 1;                   % indx will hold indices 1..K for each group
for k = 1:K                % Now we are ordering the indices for each group.
    indx = [indx; Nr(cnt:cnt+M-1)];   % So K=1 indices are in spots 1..151, K=2 in 152..303, etc
    cnt = cnt + M;
end
leftOver = Nr(cnt:end);    % We return the ordered indices and leftover values
