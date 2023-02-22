function [estimates,  model] = icatb_myquadfun(V,TR)

icatb_defaults;
global DESPIKE_OPTIONS;
global DESPIKE_SOLVER;
despike_solver = '';
if (isempty(DESPIKE_SOLVER))
    try
        despike_solver = DESPIKE_OPTIONS.solver;
    catch
    end
else
    despike_solver = DESPIKE_SOLVER;
end


if (isempty(despike_solver))
    despike_solver = 'lsqcurvefit';
end

%TR = 2;
nt = length(V);

%numP = floor(nt/30);
t = 0:TR:(nt-1)*TR;
t = t(:);
[Bhat] = icatb_regress(V,[ones(nt,1) t]);
start_point = [0 Bhat(2) Bhat(1)];
model = @myspfun;


if (strcmpi(despike_solver, 'lsqcurvefit'))
    model = @(p, t) myspfun2(p, t);
    try
        options = optimset('MaxFunEvals', 100000, 'Display', 'off', 'Algorithm','levenberg-marquardt');
        estimates = lsqcurvefit(@(p, t) myspfun2(p, t), start_point, t, V, [], [], options);
    catch
        options = {'MaxFunEvals', 100000, 'Display', 'off'};
        estimates = icatb_lsqcurvefit(@(p, t) myspfun2(p, t), start_point, t, V, [], [], options);
    end
elseif (strcmpi(despike_solver, 'fminsearch'))
    options = optimset('MaxFunEvals', 100000, 'Display', 'off');
    estimates = fminsearch(model, start_point, options);
else
    % fminunc
    options = optimset('fminunc');
    options = optimset(options, 'LargeScale', 'off', 'MaxFunEvals', 100000, 'Display', 'off');
    estimates = fminunc(@(x) myspfun(x), start_point, options);
end


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

    function yfun = myspfun2(params, t)
        % despike AFNI method; estimate a smoothish curve
        
        x0 = params;
        
        % smoothish quadratic fit
        yfun = x0(1)*t.^2 + x0(2)*t + x0(3);
        %         for ii = 1:numP
        %             yfun = yfun + x0(2+ii) * sin(2*pi*ii*t/(nt*TR)) + x0(2+ii) *cos(2*pi*ii*t/(nt*TR));
        %         end
        
        
    end

end