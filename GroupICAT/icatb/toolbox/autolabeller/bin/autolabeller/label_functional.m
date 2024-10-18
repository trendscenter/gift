% 
% 
% 
function func_pred = label_functional( sm_file, mask_file, threshold, networks, atlas_name, n, outpath )

    atlas = load_atlas( atlas_name );
    outpath = fullfile( outpath, 'nc' );

    % resample
    if ~exist( fullfile( outpath, ['r' atlas.path] ) )
        disp('resampling the atlas to input')
        sm_dat = icatb_spm_vol( sm_file );
        src_img = [sm_file ',1'];
        noisecloud_spm_coregister( src_img, atlas.path, [], outpath );
        [~,ff,xx] = fileparts(atlas.path);
        atlas.path = fullfile( outpath, ['r' ff xx] );
    end
    
    % load data
    sm_dat = icatb_spm_read_vols( icatb_spm_vol( sm_file ) );
    [vx, vy, vz, n_vols] = size( sm_dat );
    mask_dat = icatb_spm_read_vols( icatb_spm_vol( mask_file ) );
    atlas_dat = icatb_spm_read_vols( icatb_spm_vol( atlas.path ) );

    % reshape
    sm_dat = reshape( sm_dat, prod([vx, vy, vz]), [] );
    atlas_dat = reshape( atlas_dat, prod([vx, vy, vz]), [] );

    % apply the mask
    idx = find( mask_dat );
    sm_dat = sm_dat( idx, : );
    atlas_dat = atlas_dat( idx, : );

    % threshold if needed
    if threshold
        sm_dat( abs( sm_dat ) < threshold ) = 0;
    end

    disp('masking atlas')
    atlas_V_2D = convert_atlas2d( atlas.name, atlas_dat );

    disp('computing correlation')
    func_pred = cell( n_vols, 2*(n+1) );
    corrs_ = icatb_corr( sm_dat, atlas_V_2D );

    % network flags
    if networks == 1
        networks = ones( n_vols, 1 );
    end
    func_pred(:, 2) = num2cell( networks );
    
    for jj = 1:n_vols
        func_pred{jj, 1} = jj;
        [vv, ii] = maxk( corrs_(jj, :), n );
        t1 = 3;
        for kk = 1:n
            func_pred{jj, t1} = get_atlas_label( ii(kk), atlas );
            func_pred{jj, t1+1} = vv(kk);
            t1 = t1 + 2;
        end
    end
    
    % create header
    headers = {'volume', 'network'};
    for kk = 1:n
        headers = [headers, {['region_' num2str(kk)] ['spatial_corr_' num2str(kk)]}];
    end

    func_pred = [headers; func_pred];

    disp('done predicting functional labels')

function ret = load_atlas( atlas_name )
    switch atlas_name
    case 'caren'
        disp('resampling to CAREN, Doucet 2019 atlas')
        ret.name = atlas_name;
        t1 = fileparts( which( 'CAREN1_SAL.nii' ) );
        t2 = dir([t1 '/*.nii']);
        ret.path = fullfile( t2(1).folder, {t2.name} );
        ret.labels = cellfun( @(x) x(8:10), {t2.name}, 'UniformOutput', false );

    case 'gordon2016'
        disp('resampling to Gordon 2016 atlas')
        ret.name = atlas_name;
        ret.path = which( 'Parcels_MNI_111.nii' );
        labels_path = which( 'Parcels.xlsx' );
    
        t1 = readtable( labels_path );
        ret.ParcelID = t1.ParcelID;
        ret.Community = t1.Community;

    otherwise
        % yeo_buckner
        disp('resampling to Bucknerlab atlas')
        ret.name = atlas_name;
        ret.path = which( 'resampled_mask_Buckner_r286.nii' );
        labels_path = which( 'idx_286_for_Buckner17.mat' );
    
        t1 = load( labels_path );
        ret.rnames = t1.rnames;
        ret.map_ = t1.idx_286_for_Buck17.dat;        
    end

function ret = convert_atlas2d( atlas_name, atlas_dat )
    switch atlas_name
    case 'caren'
        ret = atlas_dat;
    otherwise
        s_ = [length(atlas_dat) max( atlas_dat )];
        ret = zeros(s_);
        for jj = 1:max( atlas_dat )
            idx = find( atlas_dat == jj );
            ret( idx, jj ) = 1;
        end
    end

function ret = get_atlas_label( val, atlas )
    switch atlas.name
    case 'caren'
        ret = atlas.labels{ val };
        
    case 'gordon2016'
        ret = atlas.Community{ find( atlas.ParcelID == val ) };
        
    otherwise
        % yeo_buckner
        row_ = atlas.map_( atlas.map_(:,2) == val, : );
        switch row_( 4 )
            case 1
                ret = 'subcortical';
            case 2
                ret = atlas.rnames{ row_(3) };
            case 3
                ret = 'basal ganglia';
            case 4
                ret = 'cerebellum';
        end       
    end
    
    

