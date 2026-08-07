% Only works on MAc or linux
function tests = test_icatb_backReconstruct
    tests = functiontests(localfunctions);
end

function test_backReconstruct(testCase)
    [dir_ref, fileName]=fileparts(which('gift'));
    [dir_ref, fileName]=fileparts(dir_ref);
    dirdest =[dir_ref filesep 'tests/data/test_icatb_backReconstruct'];
    cd(dirdest);
    load('sesinfo_before_br.mat');
    sesInfo.outputDir = pwd;
    sesInfo = icatb_backReconstruct(sesInfo, []);
    checkreg = load('lores_ica_br1.mat');
    match = load('lores_ica_br1_match.mat'); 
    verifyLessThan(testCase, ...
            abs(sum(sum(abs(checkreg.compSet.ic-match.compSet.ic)))), ...
            1e-10);
end

