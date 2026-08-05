function tests = test_icatb_preproc
    tests = functiontests(localfunctions);
end

function test_preproc(testCase)
    % TReNDS 8/4/26 Cyrus Eierud
    % Removes mean from each time point for fMRI
    % Verifies result is same as when running gift with default options.
    % Only works on Mac and Linux

    % Note! mat-file below contain the resulting data variable
    load('./data/test_icatb_preproc/rem_mn_timepoint.mat');

    fn = fullfile([pwd '/data/subjects/'], 'sub-007_lowres.nii');
    % Generate 50 fMRI frames as a character array
    fileN = char(arrayfun(@(k) sprintf('%s,%d', fn, k), 1:50, 'UniformOutput', false));

    preProcType='remove mean per timepoint';
    precisionType='double';
    data_orig = data;
    clear data;
    %% Preprocess data
    %
    
    % Load data
    data = icatb_read_data(fileN, [], mask_ind, precisionType);
    
    % Call pre-processing function
    modalityType = icatb_get_modality;
    if ~strcmpi(modalityType, 'conn')
        % Do not do this if connectivity matrix
        data = icatb_preproc_data(data, preProcType);
    end
    
    if (~strcmpi(preProcType, 'remove mean per timepoint'))
        % Remove mean per timepoint
        data = icatb_remove_mean(data, 1);
    end

    % Simple validation test
    verifyLessThan(testCase, sum(sum(abs(data-data_orig))), 1e-10);

end

