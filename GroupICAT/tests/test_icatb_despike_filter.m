% Only works on MAc or linux
function tests = test_icatb_despike_filter
    tests = functiontests(localfunctions);
end

function test_icatb_despikeFilter(testCase)
    [dir_ref, fileName]=fileparts(which('gift'));
    [dir_ref, fileName]=fileparts(dir_ref);
    dirdest =[dir_ref filesep 'tests/data/test_icatb_despike_filter'];
    cd(dirdest);
    load('test_despike_filter.mat');
    sesInfo.outputDir = pwd; % Set dir to user specific dir
    
    % Run the test script
    icatb_postprocess_timecourses(sesInfo, []);
    checkreg = load('nostat_postprocess_results/nostat_post_process_sub_001.mat');
    match = load('nostat_postprocess_results/nostat_post_process_sub_001_match.mat'); 
    verifyLessThan(testCase, ...
            abs(sum(sum(sum(abs(checkreg.fnc_corrs-match.fnc_corrs))))), ...
            1e-10);
    verifyLessThan(testCase, ...
            abs(sum(sum(abs(checkreg.dynamic_range-match.dynamic_range)))), ...
            1e-10);
    verifyLessThan(testCase, ...
            abs(sum(sum(abs(checkreg.fALFF-match.fALFF)))), ...
            1e-10);
    verifyLessThan(testCase, ...
            abs(sum(sum(sum(abs(checkreg.spectra_tc-match.spectra_tc))))), ...
            1e-10);
end

