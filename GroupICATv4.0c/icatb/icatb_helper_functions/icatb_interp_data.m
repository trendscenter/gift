function newTC = icatb_interp_data(tc, interpFactor)
%% Interpolate timeseries using the interpolation factor
%
% Inputs:
% 1. tc - Two dimensional matrix (T x N)
% 2. interpFactor - Interpolation factor
%
% Outputs:
% newTC - Interpolated timeseries
%

useResample = 0;
try   
    resample(1,1,2);
    useResample = 1;
catch
end

if (useResample)
    [numN, denomN] = rat(interpFactor);
    newTC = resample(tc, numN, denomN);
else
    interpFactor = ceil(interpFactor);
    for nT = 1:size(tc, 2)
        tmp2 = tc(:, nT);
        tmp2 = icatb_interp(tmp2, interpFactor);
        if (nT == 1)
            newTC = zeros(length(tmp2), size(tc, 2));
        end
        newTC(:, nT) = tmp2;
    end
    
end