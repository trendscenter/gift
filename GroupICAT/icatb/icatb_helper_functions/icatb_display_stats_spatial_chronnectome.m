function icatb_display_stats_spatial_chronnectome(schronnInfo)
% Display spatial chronnectome stats
%


load icatb_colors coldhot_sensitive;

outDir = schronnInfo.outputDir;
stats_files = schronnInfo.postprocess.stats_files;
start_terms = {};
for nR = 1:length(stats_files)
    load(fullfile(outDir, stats_files{nR}), 'spatial_trans_matrix_stats');
    UNI = spatial_trans_matrix_stats;
    for nT = 1:length(UNI.tests)
        for nC = 1:length(UNI.stats_u{nT}.Contrast)
            if (~isempty(UNI.stats_u{nT}.Contrast{nC}))
                start_terms{end + 1} = UNI.stats_u{nT}.Contrast{nC};
            else
                start_terms{end + 1} = UNI.tests{nT};
            end
        end
    end
end


[dd, inds] = unique(start_terms);

start_terms = start_terms(sort(inds));

writeV = schronnInfo.userInput.HInfo;

for ii = 1:length(stats_files)
    
    outfile = fullfile(outDir, stats_files{ii});
    load(outfile);
    [stats_dir, statsFN, extn] = fileparts(stats_files{ii});
    
    subject_clusters = cell(1, length(start_terms));
    mean_dwell_time = subject_clusters;
    spatial_trans_matrix = subject_clusters;
    
    for nTerm = 1:length(start_terms)
        
        term = start_terms{nTerm};
        
        % Subject clusters
        im_stats = zeros(prod(writeV.dim(1:3)), length(subject_cluster_stats));
        for nState = 1:length(subject_cluster_stats)
            
            UNI = subject_cluster_stats{nState};
            [stats, p, t, conname] = gatherResults(UNI, term, schronnInfo.postprocess.statsInfo);
            im_stats(schronnInfo.mask_ind, nState) = stats;
        end
        
        im_stats(~isfinite(im_stats)) = eps;
        im_stats = reshape(im_stats, [writeV.dim(1:3), length(subject_cluster_stats)]);
        out_im_name = fullfile(stats_dir, [statsFN, '_', term, '.nii']);
        icatb_write_nifti_data(fullfile(outDir, out_im_name), repmat(writeV, length(subject_cluster_stats), 1), im_stats);
        
        subject_clusters{nTerm} = out_im_name;
        
        % mean dwell time
        UNI = mean_dwell_time_stats;
        [stats, p, t, conname] = gatherResults(UNI, term, schronnInfo.postprocess.statsInfo);
        mean_dwell_time{nTerm} = t;
        
        % spatial transition matrix
        UNI = spatial_trans_matrix_stats;
        [stats, p, t, conname] = gatherResults(UNI, term, schronnInfo.postprocess.statsInfo);
        spatial_trans_matrix{nTerm} = reshape(stats, UNI.dims(2), UNI.dims(3));
        
        
    end
    
    stats_display.subject_clusters = subject_clusters;
    stats_display.mean_dwell_time = mean_dwell_time;
    stats_display.terms = start_terms;
    stats_display.spatial_trans_matrix = spatial_trans_matrix;
    
    save(outfile, 'stats_display', '-append');
    
end


%% Display results
% * *a) Subject clusters*- Statistics are done on the subject clusters for each covariate. Intensity values are reported in -sign(t)*log10(p)
% * *b) Mean dwell time* - T-values are reported for each covariate and state
% * *c) Spatial transition matrix* - Group statistics are computed on the spatial transition matrix
%
for ii = 1:length(stats_files)
    
    outfile = fullfile(outDir, stats_files{ii});
    load(outfile);
    
    for nT = 1:length(stats_display.terms)
        
        disp(['Comp ', num2str(ii), ' ', stats_display.terms{nT}]);
        
        icatb_image_viewer(fullfile(outDir, stats_display.subject_clusters{nT}), 'threshold', eps, 'image_values', 'positive and negative', 'slices_in_mm', ...
            (-40:4:72), 'structFile', fullfile(fileparts(which('groupica.m')), 'icatb_templates', 'ch2bet.nii'), 'colorbar_title', '-sgn(t)log10(p)', ...
            'convert_to_zscores', 'no');
        
        varNames = {'States', 'Tvalue'};
        mdwt = stats_display.mean_dwell_time{nT}(:);
        state_vec = (1:schronnInfo.postprocess.num_clusters(ii))';
        try
            T = table(state_vec, mdwt, 'VariableNames', varNames);
            disp(T);
        catch
            fprintf('%20s\t%20s\t%20s\n', varNames{:});
            fprintf('\n');
            for n = 1:length(state_vec)
                fprintf('%20s\t%20s\t\n', num2str(state_vec(n), '%d'), num2str(mdwt(n), '%0.3f'));
            end
            fprintf('\n');
        end
        
        %gH = figure('color', 'w');
        %sh = axes('units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
        tmp_sm = stats_display.spatial_trans_matrix{nT};
        strs = cellstr(num2str((1:size(tmp_sm, 1))'));
        icatb_plot_matrix(tmp_sm, strs, strs, 'title', ['Spatial transition matrix (', stats_display.terms{nT}, ')'], ...
            'cmap', coldhot_sensitive,  'colorbar_title', '-sgn(t)log10(p)', 'clim', [log10(eps), -log10(eps)], 'xlabel', 'Bins', ...
            'ylabel', 'Bins');
        %set(findobj(gH, 'type', 'axes'), 'Ycolor', 'k');
        %set(findobj(gH, 'type', 'axes'), 'Xcolor', 'k');
        
    end
    
    
    
end





function [p_masked, p]  = get_sig_pvalues(p, thresh, criteria)
% apply fdr correction

p_masked = thresh;
if (strcmpi(criteria, 'mafdr'))
    p = mafdr(p);
elseif(strcmpi(criteria, 'fdr'))
    p_masked = icatb_fdr(p, thresh);
end


function [term_no, con_no] = getTermAndConNum(term_name, UNI)
% Get term and contrast no
%

term_no = [];
con_no = [];
for nT = 1:length(UNI.tests)
    for nCon = 1:length(UNI.stats_u{nT}.Contrast)
        %tmp = [UNI.tests{nT}, ' ', UNI.stats{nT}.Contrast{nCon}];
        chk = (strcmpi(term_name, UNI.stats_u{nT}.Contrast{nCon}) || strcmpi(term_name, UNI.tests{nT}));
        if (chk)
            term_no = nT;
            con_no = nCon;
            return;
        end
    end
end

function [stats, p, t, conname] = gatherResults(UNI, term, statsInfo)
% gather results

stats = NaN(1, size(UNI.t_u{1}, 2));
p = stats;
t = stats;

[mIND, con_no] = getTermAndConNum(term, UNI);

if ~isempty(mIND)
    ps = UNI.p_u{mIND}(con_no, :);
    ts = UNI.t_u{mIND}(con_no, :);
    [p_masked, ps]  = get_sig_pvalues(ps, statsInfo.p_threshold, statsInfo.threshdesc);
    good_inds = ps < p_masked;
    stats(good_inds) = -log10(ps(good_inds) + eps).*sign(ts(good_inds));
    conname = UNI.stats_u{mIND}.Contrast{con_no};
    p(good_inds) = ps(good_inds);
    t(good_inds) = ts(good_inds);
end

if (~exist('conname', 'var'))
    conname = term;
end