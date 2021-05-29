function [tc_out, frames_replaced] = icatb_despike_tc(tc, TR, despike_method)
%% Despike timeseries
%
% Inputs:
% tc - Timecourses of dimensions time x components
% TR - Experimental TR in seconds
%
% Outputs:
% tc_out - Despiked timeseries
% frames_replaced - Timepoints that are replaced at each component
%

icatb_defaults;
global DESPIKE_OPTIONS;

if (~exist('despike_method', 'var'))
    
    despike_method = 1;
    try
        despike_method = DESPIKE_OPTIONS.method;
    catch
    end
    
end


if (numel(tc) == length(tc))
    tc = tc(:);
end

tc_out = zeros(size(tc));
frames_replaced = cell(1, size(tc, 2));


if (despike_method == 1)
    for ncols = 1:size(tc, 2)
        [tc_out(:, ncols), frames_replaced{ncols}] = despikeTC(tc(:, ncols), TR);
    end
elseif (despike_method == 2)
    for ncols = 1:size(tc, 2)
        [tc_out(:, ncols), frames_replaced{ncols}] = despikeTC_smooth(tc(:, ncols));
    end
else
    for ncols = 1:size(tc, 2)
        [tc_out(:, ncols), frames_replaced{ncols}] = despike_median_filter(tc(:, ncols));
    end
end


function [tc_out, ind] = despikeTC(tc, TR)

global DESPIKE_OPTIONS;


c1 = 2.5;
try
    c1 = DESPIKE_OPTIONS.c(1);
catch
end

c2 = 3;
try
    c2 = DESPIKE_OPTIONS.c(2);
catch
end
%tc = tc(:);

modelX = [ones(length(tc),1) (-1:2/(length(tc)-1):1)'];
%[lestimates] = icatb_regress(tc,[ones(length(tc),1) (-1:2/(length(tc)-1):1)']);

lestimates = (modelX'*modelX) \ (modelX'*tc);

[qestimates,  modelq] = icatb_myquadfun(tc,TR);
[splestimates,  models] = icatb_mysplinefun(tc,TR);


ylfit =  lestimates(1) + lestimates(2)*(-1:2/(length(tc)-1):1)';
yqfit = icatb_getQuadFit(qestimates,length(tc),TR);
ysfit = icatb_getSplineFit(splestimates,length(tc),TR);

err = [icatb_gfit2(tc,ylfit,'1') icatb_gfit2(tc,yqfit,'1') icatb_gfit2(tc,ysfit,'1')];

[mnerr mnID] = min(err);

if mnID == 1
    yfit =  ylfit;
elseif mnID == 2
    yfit = yqfit;
else
    yfit = ysfit;
end

res = tc - yfit;
mad_res = median(abs(res - median(res))); % median absolute deviation of residuals
sigma = mad_res* sqrt(pi/2);
s = res/sigma;
s_out = s;

ind = find(abs(s) > c1);
for uu = 1:length(ind)
    s_out(ind(uu)) = sign(s(ind(uu)))*(c1+((c2-c1)*tanh((abs(s(ind(uu)))-c1)/(c2-c1))));
end

tc_out = yfit + s_out*sigma;


function [tc_out, bad_inds] = despikeTC_smooth(tc)

global DESPIKE_OPTIONS

c1 = 2.5;
try
    c1 = DESPIKE_OPTIONS.c(1);
catch
end

c2 = 3;
try
    c2 = DESPIKE_OPTIONS.c(2);
catch
end

order = 5;
try
    order = DESPIKE_OPTIONS.order;
catch
end

if (~isempty(which('smooth.m')))
    yfit = smooth(tc, order);
else
    yfit = icatb_moving_average(tc, order);
end

res = abs(tc - yfit);
std_dev = std(res);
bad_inds = find(res > c1 * std_dev);


tc_out = tc;
if (length(bad_inds) > 0)
    tc_out(bad_inds) = yfit(bad_inds);
end

% mad_res = median(abs(res - median(res))); % median absolute deviation of residuals
% sigma = mad_res* sqrt(pi/2);
% s = res/sigma;
% s_out = s;
%
% ind = find(abs(s) > c1);
% for uu = 1:length(ind)
%     s_out(ind(uu)) = sign(s(ind(uu)))*(c1+((c2-c1)*tanh((abs(s(ind(uu)))-c1)/(c2-c1))));
% end
%
% tc_out = yfit + s_out*sigma;
%
%
% bad_inds=ind;



function [tc_out, bad_inds] = despike_median_filter(tc)

global DESPIKE_OPTIONS;

median_order = 3;
try
    median_order = DESPIKE_OPTIONS.order;
catch
end


tc_out = medfilt1(tc, median_order);

% guess frames replaced
bad_inds = find(abs(tc_out - tc) > 1e-6);