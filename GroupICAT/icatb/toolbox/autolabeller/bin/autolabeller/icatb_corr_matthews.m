function correlations = icatb_corr_matthews( label, prediction )
    % todo test

    % label = label(:);
    % perform columnwise Matthew's correlation
    
    if ~ isequal( size(label, 1), size( prediction, 1 ) ) 
        error('label does not have the same number of elements as prediction rows.')
    end

    ncol = size( label, 2 );
    nroi = size( prediction, 2 );
    correlations = NaN * ones( ncol, nroi );
    if isempty(gcp('nocreate'))
        parpool( maxNumCompThreads )
    end
    for ii = 1:ncol
        fprintf('%d ', ii)
        cm = NaN*ones(2, 2, nroi);
        parfor j = 1:nroi
            cm(:, :, j) = confusionmat( label(:, ii), prediction(:,j) );
        end
    
        tn = squeeze( cm(1, 1, :) )';
        tp = squeeze( cm(2, 2, :) )';
        fp = squeeze( cm(1, 2, :) )';
        fn = squeeze( cm(2, 1, :) )';
    
        % compute Matthews correlation coefficient (MCC)
        t1 = (tp.*tn - fp.*fn) ./ sqrt( (tp + fp) .* (tp + fn) .* (tn + fp) .* (tn + fn) );
        correlations( ii, : ) = abs( t1 );
    end
    delete(gcp('nocreate'))

end
