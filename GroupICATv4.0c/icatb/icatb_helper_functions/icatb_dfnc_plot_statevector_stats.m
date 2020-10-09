function G = icatb_dfnc_plot_statevector_stats(k, F, TM, MDT, NT, G)

if (~exist('G', 'var'))
    G = figure;
end


if isvector(F)
    SS = 1; % single subject
else
    SS = 0; % group
end

subplot(2,2,1)
if SS
    plot(1:k, F, 'm');
else
    icatb_plot_with_ste_area(gca, 1:k, F, [], 'm');
end
box off; set(gca, 'TickDir', 'out')
ylabel('Frequency')
xlabel('State (cluster index)')


subplot(2,2,3)
if SS
    plot(1:k, MDT, 'm');
else
    icatb_plot_with_ste_area(gca, 1:k, MDT, [], 'm');
end
box off; set(gca, 'TickDir', 'out')
ylabel('Mean dwell time (windows)')
xlabel('State (cluster index)')


subplot(2,2,[2,4])
if SS
    imagesc(TM, [0 1]);
else
    imagesc(squeeze(mean(TM)), [0 1])
end
colormap(jet);
C = colorbar;
set(get(C, 'YLabel'), 'String', 'Probability')
axis square
set(gca, 'XTick', 1:k, 'YTick', 1:k)
xlabel('State at t')
ylabel('State at t-1')
if SS
    title(sprintf('Number of transitions: %d', NT))
else
    title(sprintf('Number of transitions: %0.1f +/- %0.1f', mean(NT), std(NT)))
end


