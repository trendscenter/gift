% Only works on Mac or linux
% Cyrus 082426
function tests = test_icatb_nmark_batch
    tests = functiontests(localfunctions);
end

function test_icatb_nmark_batch_inner(testCase)
    [dir_ref, ~]=fileparts(which('gift'));
    [dir_ref, fileName]=fileparts(dir_ref);
    dirdest =[dir_ref filesep 'tests/data/test_icatb_nmark_batch'];
    cd(dirdest);

    % Run the test script
    icatb_batch_file_run('input_nmark_rsn.m'); 
    % after script the dir is changed to output dir (tempdir)

    % batch prefix should be nmarkrsn, yielding nmarkrsn_ica_br1.mat
    load('nmarkrsn_ica_br1.mat'); % reading struct compSet

    % Verify results
    match_tc = load([dirdest filesep 'nmarkrsn_ica_br_tc.mat']);
    match_ic = load([dirdest filesep 'nmarkrsn_ica_br_ic.mat']); 
    verifyLessThan(testCase, abs(sum(sum(abs(compSet.ic - match_ic.nmarkrsn_ica_br_ic)))), 1e-10);
    verifyLessThan(testCase, abs(sum(sum(abs(compSet.tc-match_tc.nmarkrsn_ica_br_tc)))), 1e-10);
end

