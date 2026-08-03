function tests = test_icatb_calculate_pca
    tests = functiontests(localfunctions);
end

function testDeterministicPCA(testCase)
    % Deterministic regression test for icatb_calculate_pca.
    % Verifies that PCA output remains unchanged for a fixed random seed.
    
    rng(42,'twister');
    data=rand(220,10000);
    data=data';
    pca_opts = struct('stack_data', 'yes', ...
       'storage', 'full', 'precision', 'double', ...
       'eig_solver', 'selective' );
    [pcasig, dewhiteM, Lambda, V, whiteM] = ...
       icatb_calculate_pca(data, 25, 'type', ...
       'standard', 'whiten', 1, 'verbose', 1, ...
       'preproc_type', 'none', ...
       'pca_options', pca_opts, 'remove_mean', 1);
    
    % Simple validation test
    expected_val = 1.1637705649918482;
    verifyLessThan(testCase, ...
            abs(abs(pcasig(1,1))-expected_val), ...
            1e-10);
end

