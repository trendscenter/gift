% input: path to nifti file, index of networks and correlations with how many regions
% anatomical regions to return
% output: volumes x 2*(n+1) cell of indices, anatomical region names and correlation values
function [anat_pred, corrs_] = label_anatomical( sm_file, mask_file, threshold, networks, atlas, n, outpath )
    aal_filename = 'aal2mni152.nii';
    aal_path = which( aal_filename );
    aal_label_path = which('aal.nii.txt');
    outpath = fullfile( outpath, 'nc' );

    % load SPM anatomical labels
    aal_labels = readtable( aal_label_path );

    % resample
    disp('resampling the AAL atlas to input')
    sm_dat = icatb_spm_vol( sm_file );
    src_img = [sm_file ',1'];
    noisecloud_spm_coregister( src_img, aal_path, [], outpath );

    % load and flatten all data 
    sm_dat = icatb_spm_read_vols( icatb_spm_vol( sm_file ) );
    [vx, vy, vz, n_vols] = size( sm_dat );
    mask_dat = icatb_spm_read_vols( icatb_spm_vol( mask_file ) );
    [mx, my, mz] = size( mask_dat );
    aal_dat = icatb_spm_read_vols( icatb_spm_vol( fullfile( outpath, ['r' aal_filename] ) ) );

    % check if the mask has correct dimension
    assert( isequal( [vx, vy, vz], [mx, my, mz] ), 'The mask has a different dimension than the spatial map(s).' )

    % reshape
    sm_dat = reshape( sm_dat, prod([vx, vy, vz]), [] );
    aal_dat = reshape( aal_dat, prod([vx, vy, vz]), [] );

    % apply the mask
    idx = find( mask_dat );
    sm_dat = sm_dat( idx, : );
    aal_dat = aal_dat( idx, : );

    disp('masking AAL atlas')
    s_ = [length(aal_dat) max( aal_dat )];
    aal_V_4D = zeros(s_);
    for jj = 1:max( aal_dat )
        idx = find( aal_dat == jj );
        aal_V_4D( idx, jj ) = 1;
    end

    % threshold if needed
    if threshold
        % old logic
        sm_dat( abs( sm_dat ) < threshold ) = 0;

        % % find sm_dat < threshold
        % idx_sm = find( abs( sm_dat(:) ) < threshold );
        % % set sm_dat < threshold to zero
        % % sm_dat(:) = 1;
        % sm_dat( idx_sm ) = 0;
        % % find common sm_dat < threshold and atlas = 0
        % mask_ = intersect( idx_sm, find( ~aal_V_4D(:) ) );
        % % set mask to NaN in both
        % sm_dat( mask_ ) = NaN;
        % aal_V_4D( mask_ ) = NaN;
    end

    % network flags
    anat_pred = cell( n_vols,  2*(n+1) );
    if networks == 1
        networks = ones( n_vols, 1 );
    end
    anat_pred(:, 2) = num2cell( networks );

    disp('computing correlation')
    corrs_ = icatb_corr( sm_dat, aal_V_4D );
    % % perform Matthew's correlation
    % corrs_ = icatb_corr_matthews( sm_dat, aal_V_4D );

    headers = {'volume', 'network'};
    header_flag = 0;
    for jj = 1:n_vols
        anat_pred{jj, 1} = jj;
        [vv, ii] = maxk( corrs_(jj, :), n );
        t1 = 3;
        for kk = 1:n
            anat_pred{jj, t1} = aal_labels{ aal_labels.Var1 == ii(kk), 'Var2' }{1};
            anat_pred{jj, t1+1} = vv(kk);
            t1 = t1 + 2;
            if ~header_flag
                headers = [headers, {['region_' num2str(kk)] ['spatial_corr_' num2str(kk)]}];
            end
        end
        header_flag = 1;
    end

    anat_pred = [headers; anat_pred];

    disp('done predicting anatomical labels')


