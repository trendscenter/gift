% Only works on MAc or linux
function tests = test_icatb_preproc_var
    tests = functiontests(localfunctions);
end

function test_icatb_preproc_var_inner(testCase)
    [dir_ref, fileName]=fileparts(which('gift'));
    [dir_ref, fileName]=fileparts(dir_ref);
    dirdest =[dir_ref filesep 'tests/data/test_icatb_preproc_var'];
    cd(dirdest);
    load('data_sub7_before.mat');
    
    % Run the test script
    preproc_data_remmnptp_post_result = icatb_preproc_data(data_test_pre_varnorm, 'remove mean per timepoint');
    preproc_data_varnorm_post_result = icatb_preproc_data(data_test_pre_varnorm, 'variance normalization');
    match1 = load('preproc_data_remmnptp_post.mat'); % match1.preproc_data_remmnptp_post
    match2 = load('preproc_data_varnorm_post.mat'); % match2.preproc_data_varnorm_post
    verifyLessThan(testCase, ...
            abs(sum(sum(abs(preproc_data_remmnptp_post_result-match1.preproc_data_remmnptp_post)))), ...
            1e-10);
    verifyLessThan(testCase, ...
            abs(sum(sum(abs(preproc_data_varnorm_post_result-match2.preproc_data_varnorm_post)))), ...
            1e-10);
end

