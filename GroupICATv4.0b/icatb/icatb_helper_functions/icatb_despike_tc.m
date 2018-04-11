function [tc_out, frames_replaced] = icatb_despike_tc(tc, TR)
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

if (numel(tc) == length(tc))
    tc = tc(:);
end

tc_out = zeros(size(tc));
frames_replaced = cell(1, size(tc, 2));
for ncols = 1:size(tc, 2)
    [tc_out(:, ncols), frames_replaced{ncols}] = despikeTC(tc(:, ncols), TR);
end



function [tc_out, ind] = despikeTC(tc, TR)

c1 = 2.5;
c2 = 3;
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
