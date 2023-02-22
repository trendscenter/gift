function [L,r] = icatb_test_randomness(inpMat, matType, nIter)
%% Test randomness of matrix
%

isGUI = 0;
if (~exist('inpMat', 'var') || isempty(inpMat))
    isGUI= 1;
    inpMat= icatb_selectEntry('typeEntity', 'file', 'typeSelection', 'multiple', 'filter', '*.mat', 'title', 'Select a matrix file ...');
end

drawnow;

if (isempty(inpMat))
    error('Input matrix is not selected');
end


if (ischar(inpMat))
    inpMat= icatb_load_ascii_or_mat(inpMat);
end

if (isGUI)
    numParameters = 1;
    matrix_opts = {'Non Symmetric', 'Symmetric', 'Correlation'};
    inputText(numParameters).promptString = 'Select Matrix Type';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = matrix_opts;
    inputText(numParameters).dataType = 'string';
    inputText(numParameters).tag = 'mat_type';
    inputText(numParameters).enable = 'on';
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = 'Enter number of iterations to estimate means and covariance of singular values';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = num2str(10000);
    inputText(numParameters).dataType = 'numeric';
    inputText(numParameters).tag = 'max_iter';
    inputText(numParameters).enable = 'on';
    
    answers = icatb_inputDialog('inputtext', inputText, 'Title', 'Enter parameters for estimating randomness', 'handle_visibility',  'on', 'windowStyle', 'modal');
    
    chk = strcmpi(matrix_opts, answers{1});
    matType = double(find(chk == 1)) - 1;
    
    nIter = answers{2};    
    
end

drawnow;

disp('Computing randomness metric L ...');
[L, r] = icatb_dist2RandMat(inpMat, matType, nIter);
disp('Done');

H=figure;
imagesc(inpMat);
title(['L=', sprintf('%.0f',r.L), '  pval=',sprintf('%.2e',r.pval)]);
