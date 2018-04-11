function [estimates,  model] = icatb_myquadfun(V,TR)
%TR = 2;
nt = length(V);

%numP = floor(nt/30);
t = 0:TR:(nt-1)*TR;
t = t(:);
[Bhat] = icatb_regress(V,[ones(nt,1) t]);
start_point = [0 Bhat(2) Bhat(1)];
model = @myspfun;
options = optimset('MaxFunEvals', 100000, 'Display', 'off');
estimates = fminsearch(model, start_point, options);

    function [err yfun] = myspfun(params)
        % despike AFNI method; estimate a smoothish curve

        x0 = params;

        % smoothish quadratic fit
        yfun = x0(1)*t.^2 + x0(2)*t + x0(3);
        %         for ii = 1:numP
        %             yfun = yfun + x0(2+ii) * sin(2*pi*ii*t/(nt*TR)) + x0(2+ii) *cos(2*pi*ii*t/(nt*TR));
        %         end

        err = sum((V - yfun).^2);

    end

end