% noisecloud_download is intended for only downloading feature scripts from
% the noisecloud database into the temp folder

function noisecloud_download

    %% Step 1: GET FEATURES FROM noisecloud for download
    fprintf('%s\n','Reading spatial and temporal features from noisecloud...');
    spatial = urlread('http://vbmis.com/bmi/ncdb/rest/spatial');
    temporal = urlread('http://vbmis.com/bmi/ncdb/rest/temporal');

    % Parse spatial and temporal features
    fprintf('%s\n','Parsing feature data...');
    % Parse json for subjects and concepts
    [ spatial ~ ] = noisecloud_parse_json(spatial); spatial = spatial{1};
    [ temporal ~ ] = noisecloud_parse_json(temporal); temporal = temporal{1};

    feature_labels = [];
    idx = 1;

    % Write temporary scripts for spatial and temporal features
    for s=1:size(spatial,2)
        nfeatures = str2num(spatial{s}.n);
        labels = regexp(spatial{s}.label,',','split');
        if nfeatures ~= size(labels,2) % If there is one label to describe features
            for f=1:nfeatures
                feature_labels{idx} = [ spatial{s}.label '.' num2str(f) ];
                idx = idx +1;
            end
        else % If there is a unique label per feature
            for f=1:nfeatures
                feature_labels{idx} = labels{f};
                idx = idx +1;
            end
        end
        script = fopen([ 'temp/' spatial{s}.name '.m' ],'w');
        lines = regexp(spatial{s}.code,';','split');
        fprintf(script,'%s',lines{1}); 
        for l=2:size(lines,2)
           fprintf(script,'%s\n%s',';',[ lines{l} ]); 
        end
        fclose(script);
    end

    for t=1:size(temporal,2)
        nfeatures = str2num(temporal{t}.n);
        labels = regexp(temporal{t}.label,',','split');
        if nfeatures ~= size(labels,2) % If there is one label to describe features
            for f=1:nfeatures
                feature_labels{idx} = [ temporal{t}.label '.' num2str(f) ];
                idx = idx + 1;
            end
        else % If there is a unique label per feature
            for f=1:nfeatures
                    feature_labels{idx} = labels{f};
                idx = idx + 1;
            end
        end
        script = fopen([ 'temp/' temporal{t}.name '.m' ],'w');
        lines = regexp(spatial{s}.code,';','split');
        fprintf(script,'%s',lines{1}); 
        for l=2:size(lines,2)
            fprintf(script,'%s\n%s',';',[ lines{l} ]); 
        end
        fclose(script);
    end

    fprintf('%s\n','Done downloading features.  Scripts can be found in the temp folder.');

end