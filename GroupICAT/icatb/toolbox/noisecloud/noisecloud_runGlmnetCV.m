% noisecloud_runGlmnetCV() will perform K fold cross validation with logistic 
% regression on a set of input data, and produce a cross validated
% ROC curve, and a P value from bootstrapping a user specified number of
% times.
%
% USAGE: runglmnetCV(X,y,K,alphav,Lambda,type,iterations)
% INPUT: 
%        X: an n x m matrix of n objects (rows) each with m features (cols)
%        y: binary labels (0,1) to be changed to (1,2) 
%        K: the number of cross fold validations
%        iterations: to use for bootstrapping (to get p value)
%        alphav: a range of alpha values to do grid search to find optimum
%        Lambda: try starting with a range, then just using optimal, eg:
%                [0.0625 0.09 0.125 0.25 0.5 1 2 4 8 16 32]
%        type: 'binomial' only supported / tested
%        colly_labels: correspond to feature labels
%
% SUGGESTED:      K =10;
%                 alphav = 1:-0.01:0.1;
%                 Lambda = 0.0625;
%                 type = 'binomial';

function REALRESULT = noisecloud_runGlmnetCV(X,y,K,iterations,alphav,Lambda,type,colly_labels)

    % Run it once to get a RESULT
    REALRESULT = noisecloud_glmnetCV(X,y,K,Lambda,alphav,type);

    % Figure out which features are weighted most highly
    % NOTE that REALRESULT.wopt is +1 in length LONGER than the number of features BECAUSE
    % the first weight is the intercept value (see glmnetPredict nbeta
    % Find highest weighted features
    [weights,indexy] = sort(REALRESULT.wopt(2:end),1,'descend');

    % Then find the indices where the sorted weights aren't zero
    REALRESULT.nonzeroFeatures = colly_labels(indexy(weights ~= 0));
    REALRESULT.nonzeroWeights = weights(weights ~= 0);
    
    acc_results = zeros(iterations,1);

    if iterations > 1
            fprintf('\n%s\n',[ 'Running ' num2str(iterations) ' permutations with randomly shuffled labels to get a p value' ]);
            f = figure; subplot(1,2,1); hold on
            
            % Plot red line for random chance
            line(0:1,0:1,'color','r')
            title([ 'ROC Curve for ' num2str(iterations) ' Permutations' ]);
            for it=1:iterations
                % Create random labels
                yrand = y(randperm(length(y)));
                PERMRESULT = noisecloud_glmnetCVperm(X,yrand,K,Lambda,alphav,type);
                acc_results(it) = PERMRESULT.best_cva;
                
                % Add sensitivity and specificity to plot
                plot(PERMRESULT.spec,PERMRESULT.sens,'LineWidth',2);
                
                % Confusion Matrix  [TP FP; FN TN]
                subplot(1,2,2);
                data = [PERMRESULT.EVAL(8) PERMRESULT.EVAL(10); PERMRESULT.EVAL(11) PERMRESULT.EVAL(9) ];
                cnames = {'YES','NO'};
                rnames = {'YES','NO'};
                set(gca,'Units','normalized');
                t = uitable('Data',data,'ColumnName',cnames,'RowName',rnames,'Tag','Table 4','TooltipString','Table 4','Parent',f,'Units','Normalized','Position',[.5 .5 .5 .3]);
                title('Confusion Matrix'); axis off
                subplot(1,2,1); title([ 'Permutation ' num2str(it) ' : Best CVA ' num2str(PERMRESULT.best_cva) ]);
            end

    % Normalize our distribution
    acc_std = std(acc_results);
    acc_mean = mean(acc_results);
    % acc_results_norm = (acc_results - acc_mean) / acc_std;

    % If you want to do statistical calculations for the permutations, add
    % them here!
    % p value and confidence interval for 95% of distribution
    % [~,p,ci] = ttest(acc_results,5);
    REALRESULT.perm.mean = acc_mean;
    REALRESULT.perm.n = iterations;
    REALRESULT.perm.std = acc_std;

    end
    

    % CONFIDENCE INTERVALS
    % We can calculate 95% confidence intervals for our results based on
    % the standard error of a binomial.  Given k successes out of n trials, 
    % the observered proportion p = k/n, with standard error SE = sqrt(p*(1-p)/n).  
    % The 95% confidence interval is p +/- 1.96*SE.  
    % If n < 20, you may need to use an adjusted formula.

    % CVA confidence interval
    % P is our best cross validation accuracy, k/n
    p = REALRESULT.best_cva;
    SE = sqrt(p*(1-p)/size(X,1));
    int95lower = p - 1.96*SE;
    int95upper = p + 1.96*SE;
    REALRESULT.CI95.best_cva = [int95lower int95upper];
    
    % Sensitivity confidence interval
    kk = REALRESULT.ROC.tp;
    n = sum(y == 1);  % noise is 1, not is 2 
    p = kk/n;
    SE = sqrt(p*(1-p)/n);
    int95lower = p - 1.96*SE;
    int95upper = p + 1.96*SE;
    REALRESULT.CI95.sens = [int95lower int95upper];
    
    % Specificity confidence interval
    kk = REALRESULT.ROC.tn;
    n = sum(y == 2);  % noise is 1, not is 2 
    p = kk/n;
    SE = sqrt(p*(1-p)/n);
    int95lower = p - 1.96*SE;
    int95upper = p + 1.96*SE;
    REALRESULT.CI95.spec = [int95lower int95upper];

end
