function icatb_plot_dyn_coh_clusters(Centroids, idx, coin, comps, TR)

k_cluster = size(Centroids, 1);

incoi_mask_v = find(coin);
incoi_mask = coin;

numComps = length(comps);

alaki_mask = ones(numComps);
mask_inds = alaki_mask - triu(alaki_mask);


idx = reshape(idx, [], length(coin(:)));
[ time_freq_inds_unwarped ] = UnwarpWavelet(idx(:, incoi_mask_v), incoi_mask);
%time_freq_inds_unwarped = best_cluster_results;

bincounts = histc(time_freq_inds_unwarped(:),min(time_freq_inds_unwarped(:)):max(time_freq_inds_unwarped(:)));
%bincounts = bincounts(2:end);
[sort_vals,sort_inds] = sort(bincounts,'descend');

% Centroids = best_cluster_centroids;

Masked_True_Centroids = Centroids;

True_Centroids = zeros(size(Centroids,1), prod(size(mask_inds)));
True_Centroids(:,find(mask_inds)) = Centroids;
True_Centroids = reshape(True_Centroids, size(Centroids,1), size(mask_inds,1), size(mask_inds,2));
True_Centroids(isnan(True_Centroids))=0;


sz = get(0, 'screensize');

figW = round(sz(3)*0.82);
figH = round(sz(4)*0.82);

figPos = [0.5*sz(3) - 0.5*figW, 0.5*sz(4)  - 0.5*figH, figW, figH];


for curr_cmp_ind = 1:k_cluster
    
    F= figure('color', 'w', 'name', ['State ', num2str(curr_cmp_ind)], 'tag', ['state_dyn_', num2str(curr_cmp_ind)]);
    set(F, 'resize', 'on');
    set(F, 'position', figPos);
    
    sh = subplot(2, 3, [1,2,4,5]);
    
    graphicsH(curr_cmp_ind).H = F;
    
    
    curr_cmp = squeeze(True_Centroids(sort_inds(curr_cmp_ind),:,:));
    curr_cmp = curr_cmp + rot90(flipud(curr_cmp),-1);
    curr_centroid_inds = idx == sort_inds(curr_cmp_ind);
    subject_specific_occupancy = sum(curr_centroid_inds,2);
    msgStr = sprintf('State %d appears in %1.2f %% of Subjects', curr_cmp_ind, 100*length(find(subject_specific_occupancy > 0))/length(subject_specific_occupancy));
    
    fprintf(msgStr);
    fprintf('\n');
    
    %subplot(max_clusters,max_clusters-start_cluster+1,(max_clusters-start_cluster+1)*(curr_cmp_ind-1)+k_cluster-start_cluster+1);
    %subplot(k_cluster,1,curr_cmp_ind);
    %     cbfreeze;
    %     freezeColors;
    %     subplot(2,2,1);
    
    
    %     load('MyColormaps_b','mycmap')
    %     CM_over = mycmap;
    
    CM_over = hsv(30);%colormap('hsv');
    CM_over = circshift(CM_over,round(length(CM_over)*0.5));
    icatb_magphase_plot_labels(curr_cmp,CM_over);axis ij%colorbar;
    title(msgStr);
    set(sh, 'XTickLabel', num2str(comps(:)));
    set(sh, 'YTickLabel', num2str(comps(:)));
    xlabel('Components');
    ylabel('Components');
    
    
    sh = subplot(2, 3, 3);
    
    %figure;
    time_freq_inds_tmp = (time_freq_inds_unwarped == sort_inds(curr_cmp_ind));
    time_freq_inds_tmp = permute(time_freq_inds_tmp, [2 1 3]);
    time_freq_inds_tmp = reshape(time_freq_inds_tmp, size(time_freq_inds_tmp,1), []);
    time_freq_inds_tmp = sum(time_freq_inds_tmp, 2);
    if (curr_cmp_ind == 1)
        freq_hist = linspace(0.01, 1/TR/2, length(time_freq_inds_tmp));
        freq_hist_str = cellstr(num2str(freq_hist(:),'%0.2f'));
    end
    bar(time_freq_inds_tmp);
    axis(sh, 'square');
    set(sh, 'xticklabel', freq_hist_str);
    title('Frequency Histogram', 'parent', sh);
    xlabel('Hz', 'parent', sh);
    
    
    sh = subplot(2, 3, 6);
    plotPhaseHist(CM_over, sh);
    
    
end


icatb_plotNextPreviousExitButtons(graphicsH);



function [ UnWarped_Output ] = UnwarpWavelet( Warped_Input, Mask )
%UNWARPWAVELET Summary of this function goes here
%   Detailed explanation goes here

incoi_mask_v = find(Mask > 0);


UnWarped_Output = -ones(size(Warped_Input,1), size(Mask,1), size(Mask,2));


for curr_step = 1:size(Warped_Input,1)
    tmpUnWarped = zeros(size(Mask));
    tmpUnWarped(incoi_mask_v) = Warped_Input(curr_step,:);
    
    for curr_row = 1:size(Mask,1)
        inXInds = find(Mask(curr_row,:) == 1);
        if(~isempty(inXInds))
            inX = ((size(Mask,2)-1)/(inXInds(end)-inXInds(1)))*(inXInds - inXInds(1)) + 1;
            inY = tmpUnWarped(curr_row, inXInds);
            UnWarped_Output(curr_step, curr_row, :) = interp1(inX(:), inY(:), 1:size(Mask,2), 'nearest','extrap');
        end
    end
    
end


function plotPhaseHist(cmap, sh)
% Plot phase histogram

startAngle = 0;
numSlices = size(cmap, 1);
pi_incr = 2*pi/numSlices;
rin = 1;
rout = 1.3;
rText = 1.8;

chkP = ceil(linspace(0, numSlices, 5));
textsToPlot = {'0', '+ \pi /2', '+/- \pi', '-\pi /2'};
ths = linspace(0, 2*pi, length(chkP));

for n = 1:numSlices
    %tmpC = colorFrames{n};
    endAngle = startAngle + pi_incr;
    t2 = linspace(startAngle, endAngle, 256);
    startAngle = endAngle;
    xin = rin*cos(t2(end:-1:1));
    yin = rin*sin(t2(end:-1:1));
    xout = rout*cos(t2);
    yout = rout*sin(t2);
    
    patch([xout, xin],[yout, yin], cmap(n,:), 'edgecolor', 'none', 'parent', sh);
    hold(sh, 'on');
end

axis(sh, 'equal');

for n = 1:length(textsToPlot)
    
    xt = rText*cos(ths(n));
    yt = rText*sin(ths(n));
    text(xt, yt, textsToPlot{n}, 'fontsize', 11, 'fontweight', 'bold', 'interpreter', 'tex');
    
end

set(sh, 'box', 'off');
set(sh, 'xtick', []);
set(sh, 'ytick', []);
set(sh,'xcolor', get(get(sh, 'parent'), 'color'));
set(sh,'ycolor', get(get(sh, 'parent'), 'color'));
