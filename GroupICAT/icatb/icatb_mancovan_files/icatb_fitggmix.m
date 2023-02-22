function [y, cutoff, PF, P0, FH] = icatb_fitggmix(meanmap, nstd)
%[y, cutoff, PF, P0, FH] = fitggmix(meanmap, nstd)
[N,X] = hist(meanmap,800);


XDATA = X;
YDATA = N;
area = trapz(XDATA,YDATA);
YDATA = YDATA/area;

options = {'TolFun',1e-14, 'TolX', 1e-14, 'MaxFunEvals', 100000, 'Display', 'off'};

sigma = std(meanmap);

beta = std(meanmap)*.4;
alpha = 5;
[mv, ind]=max(YDATA);
offset = XDATA(ind);

P0 =    [0.8,                 sigma,         .1,                 alpha,           beta,       alpha      beta  offset];
LB = [0.5                      0          0                 0             0          0                  0     -0.5*sigma];
UB = [1                       Inf         0.5             Inf           Inf         Inf                Inf     0.5*sigma];

try
    options = optimset(options{:});
    [PF,RESNORM,RESIDUAL,EXITFLAG,OUTPUT] = lsqcurvefit(@icatb_ggmix,P0,XDATA,YDATA, LB, UB, options);
catch
    options = {'TolFun',1e-14, 'TolX', 1e-14, 'MaxFunEvals', 100000, 'Display', 'off', 'Jacobian', 'off'};
    [PF,RESNORM,RESIDUAL,EXITFLAG,OUTPUT] = icatb_lsqcurvefit(@icatb_ggmix,P0,XDATA,YDATA, LB, UB, options);
end

%[PF,RESNORM,RESIDUAL,EXITFLAG,OUTPUT] = lsqcurvefit(@ggmix,P0,XDATA,YDATA, [],[], options);

%% plotting
%x = linspace(-max(abs(XDATA)), max(abs(XDATA)), 1000);
x = linspace(-8*PF(2), 8*PF(2), 2000);
[z, y]=icatb_ggmix(PF,x);
FH = figure; set(FH, 'Color', [0 0 0])
plot(XDATA,YDATA, 'Color', [.8 .8 .8], 'LineWidth', 1);
hold on;
plot(x,y(1,:), 'c',  x,y(2,:), 'g',  x,y(3,:), 'y')
plot(x,z,'r')
L = legend('data', 'gaussian', 'gamma_+', 'gamma_-', 'fit');
set(L, 'Color', [0 0 0], 'TextColor', [1 1 1])
axis tight;
ax = axis;
axis([min(x) max(x) 0 ax(4)])
H = gca;
set(H, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1]);
set(get(H, 'XLabel'), 'String','voxel value', 'FontName', 'Arial', 'FontSize', 12);
set(get(H, 'YLabel'), 'String','frequency', 'FontName', 'Arial', 'FontSize', 12);
set(H, 'FontName', 'Arial', 'FontSize', 10);
titstring = sprintf('cutoff: %d std, mu = %0.3f, sigma = %0.3f', nstd, PF(8), PF(2));
set(get(H, 'Title'), 'String', titstring, 'Color', [1 1 1],'FontName', 'Arial', 'FontSize', 12);

%axis([-2 2 0 ax(4)])

cutoff(1) = PF(8)-nstd*PF(2);
cutoff(2) = PF(8)+nstd*PF(2);

hold on
plot([cutoff(1) cutoff(1)], [0, max(YDATA)], 'm--')
plot([cutoff(2) cutoff(2)], [0, max(YDATA)], 'm--')

PF = param_to_struct(PF);
P0 = param_to_struct(P0);


function S = param_to_struct(P)
S.gauss_a = P(1);
S.gauss_mu = P(8);
S.gauss_sigma = P(2);
S.gamma1_a = P(3);
S.gamma1_alpha = P(4);
S.gamma1_beta = P(5);
S.gamma2_a = 1-P(1)-P(3);
S.gamma2_alpha = P(6);
S.gamma2_beta = P(7);