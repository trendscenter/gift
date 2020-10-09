function [ t, p, stats] = icatb_nan_mT(data, groups, covariates, term, options, missing_subs_cutoff)
%% Nan version of mT

icatb_defaults;
global MANCOVA_DEFAULTS;

if (~exist('missing_subs_cutoff', 'var'))
    missing_subs_cutoff = [];
    try
        missing_subs_cutoff = MANCOVA_DEFAULTS.MISSING_SUBJECTS_CUTOFF;
    catch
    end
    if (isempty(missing_subs_cutoff))
        missing_subs_cutoff = 90;
    end
end

missing_subs_cutoff = ceil((missing_subs_cutoff*size(data, 1))/100);

nan_cols = any(isnan(data));

good_cols = find(nan_cols == 0);
nan_cols = find(nan_cols == 1);

t = NaN(1, size(data, 2));
p = t;

if (~isempty(covariates))
    stats.Terms = covariates;
else
    stats.Terms = {term};
end

stats.X = groups;

%stats.B = NaN(length(stats.Terms), size(data, 2));

stats.SSE = NaN(1, size(data, 2));
stats.DFE = NaN(1, size(data, 2));
stats.MSE = NaN(1, size(data, 2));
stats.Term = term;



if (~isempty(good_cols))
    [t(good_cols), p(good_cols), stats1] = mT(data(:, good_cols), groups, covariates, term, options);
    stats.B = NaN(size(stats1.B, 1), size(data, 2));
    stats.B(:, good_cols) = stats1.B;
    stats.SSE(1, good_cols) = stats1.SSE;
    stats.DFE(1, good_cols) = stats1.DFE;
    stats.MSE(1, good_cols) = stats1.MSE;
    stats.Levels = stats1.Levels;
    clear stats1;
end

if (~isempty(nan_cols))
    countB = 0;
    for nC = 1:length(nan_cols)
        y = data(:, nan_cols(nC));
        chk = isnan(y);
        if (exist('missing_subs_cutoff', 'var'))
            if (length(find(chk == 1)) >= missing_subs_cutoff)
                continue;
            end
        end
        y(chk) = [];
        grp = groups;
        grp(chk, :) = [];
        
        if (length(y) == 1)
            continue;
        end
        
        if (~isempty(y))
            
            try
                [tmp_t, tmp_p, stats2] = mT(y, grp, covariates, term, {});
                if (~isreal(tmp_t))
                    continue;
                end
                countB = countB + 1;
                t(nan_cols(nC)) = tmp_t;
                p(nan_cols(nC)) = tmp_p;
                if (countB == 1)
                    if (~isfield(stats, 'B'))
                        stats.B = NaN(size(stats2.B, 1), size(data, 2));
                    end
                end
                stats.B(:, nan_cols(nC)) = stats2.B;
                stats.SSE(1, nan_cols(nC)) = stats2.SSE;
                stats.DFE(1, nan_cols(nC)) = stats2.DFE;
                stats.MSE(1, nan_cols(nC)) = stats2.MSE;
                stats.Levels = stats2.Levels;
            catch
                
            end
            
            clear stats2;
            
        end
    end
end

if (~isfield(stats, 'Levels'))
    stats.Levels = [1, 0];
end
