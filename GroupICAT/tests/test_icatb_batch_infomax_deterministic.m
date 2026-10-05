% Only works on Mac or linux
% Cyrus 082426
% example: results = runtests('test_icatb_batch_infomax_deterministic')
function tests = test_icatb_nmark_batch
    tests = functiontests(localfunctions);
end

function test_icatb_nmark_batch_inner(testCase)


    % Add a disp('WARNING, this icatb_defaults.m is manipulated by test_icatb_batch_infomax_deterministic.m') that the

    [dir_ref_root, ~]=fileparts(which('gift'));
    [dir_ref, fileName]=fileparts(dir_ref_root);
    dirdest =[dir_ref filesep 'tests/data/test_icatb_batch_infomax'];
    cd(dirdest);
    
    defaultsFile = [dir_ref_root filesep 'icatb_defaults.m'];
    
    % Create unique timestamped backup name
    timestamp = datestr(now, 'yyyymmdd_HHMMSS_FFF');
    backupFile = fullfile(dir_ref_root, ...
        ['icatb_defaults_backup_', timestamp, '_warning.m']);
    
    % Move original file to backup
    movefile(defaultsFile, backupFile);
    
    % Copy backup back to original filename
    copyfile(backupFile, defaultsFile);
    
    % Read new copy
    txt = fileread(defaultsFile);
    
    % Make modifications
    old1 = 'RAND_SHUFFLE = 1;';
    new1 = 'RAND_SHUFFLE = 0; %WARNING MODIFIED BY TEST SCRIPT';
    
    old2 = 'NORAND_DETERMINISTIC = 0;';
    new2 = 'NORAND_DETERMINISTIC = 1; %WARNING MODIFIED BY TEST SCRIPT';
    
    % Make sure expected lines exist
    assert(contains(txt, old1), ...
        'Could not find "%s" in %s', old1, defaultsFile);
    
    assert(contains(txt, old2), ...
        'Could not find "%s" in %s', old2, defaultsFile);
    
    % Replace them
    txt = strrep(txt, old1, new1);
    txt = strrep(txt, old2, new2);
    
    % Write modified file
    fid = fopen(defaultsFile, 'w');
    assert(fid ~= -1, 'Could not open %s for writing.', defaultsFile);
    
    cleanupObj = onCleanup(@() fclose(fid));
    fprintf(fid, '%s', txt);
    clear cleanupObj
    
    fprintf('Original saved as:\n  %s\n', backupFile);
    fprintf('Modified defaults file:\n  %s\n', defaultsFile);

    % Run the test script
    icatb_batch_file_run('input_infomax.m'); 
    % after script the dir is changed to output dir (tempdir)

    % Restore original defaults file
    delete(defaultsFile)
    movefile(backupFile, defaultsFile);    

    % batch prefix should be nmarkrsn, yielding nmarkrsn_ica_br1.mat
    load('infomaxrsn_ica_br1.mat'); % reading struct compSet

    % Verify results
    match = load([dirdest filesep 'ref_infomaxrsn_ica_br1.mat']);
    verifyLessThan(testCase, abs(sum(sum(abs(compSet.ic - match.compSet.ic)))), 1e-10);
    verifyLessThan(testCase, abs(sum(sum(abs(compSet.tc - match.compSet.tc)))), 1e-10);
end

