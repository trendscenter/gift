function ICA_Options = icatb_icaOptions(dataSize, algorithm_index, handle_visibility)
%% Return arguments depending on the ICA algorithm
% icatb_icaOptions uses icatb_inputDialog function to return ICA options
% ICA Options is a cell arary with strings followed by values
%
% Inputs:
% 1. dataSize - Number of components by frames
% 2. algorithm_index - ICA Algorithm number
% 3. handle_visibility - Options are 'on' or 'off'
%
% Outputs:
% ICA_Options - ICA Options

ICA_Options = {};

if ~exist('handle_visibility', 'var')
    handle_visibility = 'on';
end

if ispc
    windowStyle = 'modal';
else
    windowStyle = 'normal';
end

%% Determine the data size
numComps = dataSize(1);
frames = dataSize(2);

%% Available ICA algorithms
icaAlgo = icatb_icaAlgorithm;

if ischar(algorithm_index)
    algorithm_index = strmatch(lower(algorithm_index), lower(icaAlgo), 'exact');
    if isempty(algorithm_index)
        error('Check the algorithm name you specified');
    end
end

%% Check the number of options
if (algorithm_index <= size(icaAlgo, 1))
    % selected ica algorithm
    ica_algorithm = deblank(icaAlgo(algorithm_index , :));
else
    disp(['Presently there are ', num2str(size(icaAlgo, 1)), ' algorithms']);
    ICA_Options = {};
    return;
end


%% ICA algorithm options
switch (lower(ica_algorithm))
    
    case {'infomax', 'sdd ica'}
        %% Infomax algorithm or SDD ICA algorithm
        
        % Defaults for the inofmax algorithm
        MAX_WEIGHT = 1e8;       % guess that weights larger than this have blown up
        DEFAULT_BLOCK = floor(sqrt(frames/3));
        DEFAULT_weight = 0.0;
        DEFAULT_LRATE = 0.015/log(numComps);
        DEFAULT_ANNEALSTEP   = 0.9;
        DEFAULT_ANNEALDEG = 60;
        DEFAULT_STOP = 0.000001;
        DEFAULT_MAXSTEPS = 512;
        DEFAULT_MOMENTUM = 0.0;
        DEFAULT_EXTENDED = 0.0;
        
        if strcmpi(ica_algorithm, 'infomax')
            % dialog Title
            dlg_title = 'Select the Options for the Infomax algorithm';
        else
            % dialog Title
            dlg_title = 'Select the Options for the SDD ICA algorithm';
        end
        
        numParameters = 1;
        
        inputText(numParameters).promptString = [['Select block less than ', num2str(frames)]  ...
            [' where Default = ', num2str(DEFAULT_BLOCK)]];
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(DEFAULT_BLOCK);
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'block';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = ['Select stop where Default =', num2str(DEFAULT_STOP)];
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(DEFAULT_STOP);
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'stop';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = ['Select weight where Maxweight =', num2str(MAX_WEIGHT)];
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(DEFAULT_weight);
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'weights';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString ='Select lrate where min = 0.000001 and max = 0.1';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(DEFAULT_LRATE);
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'lrate';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = ['Select maxsteps where Default =' num2str(DEFAULT_MAXSTEPS)];
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(DEFAULT_MAXSTEPS);
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'maxsteps';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select anneal between (0 1]';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(DEFAULT_ANNEALSTEP);
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'anneal';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select annealdeg between [0 180]';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(DEFAULT_ANNEALDEG);;
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'annealdeg';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select momentum between [0 1]';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(DEFAULT_MOMENTUM);
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'momentum';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = ['Select Extended where default = ' num2str(DEFAULT_EXTENDED)];
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'0', '1'};
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'extended';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Number of Components';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = num2str(numComps);
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'ncomps';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = ' Select Posact';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'off', 'on'};
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'posact';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = ' Select Sphering';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'on', 'off'};
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'sphering';
        inputText(numParameters).enable = 'on';
        
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select Bias';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'on', 'off'};
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'bias';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select Verbose';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'on', 'off'};
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'verbose';
        inputText(numParameters).enable = 'on';
        
        
    case 'fast ica'
        %% Options for Fast ICA
        
        % Defaults for the FastICA
        epsilon = 0.0001;
        maxNumIterations  = 1000;
        maxFinetune       = 5;
        sampleSize        = 1;
        
        % dialog Title
        dlg_title = 'Select the Options for the Fast ICA algorithm';
        
        numParameters = 1;
        
        inputText(numParameters).promptString = ['Select epsilon where default = ', num2str(epsilon)];
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(epsilon);
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'epsilon';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString =  ['Select maximum iterations where default = ', num2str(maxNumIterations)];
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(maxNumIterations);
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'maxNumIterations';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = ['Select maximum finetune where default = ', num2str(maxFinetune)];
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(maxFinetune);
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'maxFinetune';
        inputText(numParameters).enable = 'on';
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = ['Select sampleSize between [0, 1] where default = ', num2str(sampleSize)];
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(sampleSize);
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'sampleSize';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Number of components';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = num2str(numComps);
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'numOfIC';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select verbose as on/off';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'on', 'off'};
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'verbose';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select approach as defl/symm';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'defl', 'symm'};
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'approach';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select finetune from one of the following:  [off, pow3, tanh, gauss, skew]';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'off', 'pow3', 'tanh', 'gauss', 'skew'};
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'finetune';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select stabilization as on/off';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'off', 'on'};
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'stabilization';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select g from one of the following: [pow3, tanh, gauss, skew]';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'pow3', 'tanh', 'gauss', 'skew'};
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'g';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select Only from the following:';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'all', 'pca', 'white'};
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'only';
        inputText(numParameters).enable = 'on';
        
    case 'semi-blind infomax'
        %% Semi-blind Infomax
        
        % Defaults for the inofmax algorithm
        MAX_WEIGHT = 1e8;       % guess that weights larger than this have blown up
        DEFAULT_BLOCK = floor(sqrt(frames/3));
        DEFAULT_weight = 0.0;
        DEFAULT_LRATE = 0.015/log(numComps);
        DEFAULT_ANNEALSTEP   = 0.9;
        DEFAULT_ANNEALDEG = 60;
        DEFAULT_STOP = 0.000001;
        DEFAULT_MAXSTEPS = 512;
        DEFAULT_MOMENTUM = 0.0;
        DEFAULT_EXTENDED = 0.0;
        
        % dialog Title
        dlg_title = 'Select the Options for the SBICA algorithm';
        
        numParameters = 1;
        
        inputText(numParameters).promptString = [['Select block less than ', num2str(frames)]  ...
            [' where Default = ', num2str(DEFAULT_BLOCK)]];
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(DEFAULT_BLOCK);
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'block';
        inputText(numParameters).enable = 'on';
        
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select correlation correction factor';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = '0.5';
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'prefs';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = ['Select stop where Default =', num2str(DEFAULT_STOP)];
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(DEFAULT_STOP);
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'stop';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select anneal between (0 1]';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(DEFAULT_ANNEALSTEP);
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'anneal';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select annealdeg between [0 180]';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(DEFAULT_ANNEALDEG);
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'annealdeg';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = ' Select Posact';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'off', 'on'};
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'posact';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select Bias';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'on', 'off'};
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'bias';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = ['Select Sphering'];
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'on', 'off'};
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'sphering';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = ['Select Extended where default = ' num2str(DEFAULT_EXTENDED)];
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'0', '1'};
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'extended';
        inputText(numParameters).enable = 'on';
        
    case 'constrained ica (spatial)'
        %% Spatial constrained ICA
        
        % dialog Title
        dlg_title = 'Select the options for the Multi-Fixed ICA algorithm';
        
        numParameters = 1;
        
        inputText(numParameters).promptString = 'Select number of loops';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = '10';
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'loopnum';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select closeness measure where maximum value is 4.8';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = '0.08';
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'threshold';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select positive constant (rho)';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = '1.0e-2';
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'rho';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select Lagrange Multiplier (mu) between [0.8 10.8]';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = '10.8';
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'mu';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select penalty parameter (gamma) where maximum value is 0.02';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = '0.02';
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'gamma';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select stopping criteria';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = '0.0001';
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'epsilon';
        inputText(numParameters).enable = 'on';
        
    case {'fbss', 'erbm'}
        %% FBSS
        
        dlg_title = 'Select the options for the FBSS algorithm';
        
        numParameters = 1;
        
        inputText(numParameters).promptString = 'Enter filter length';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(min(11, frames/50));
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'filter_length';
        inputText(numParameters).enable = 'on';
        
    case 'iva-gl'
        %% IVA
        
        dlg_title = 'Select the options for the IVA-GL algorithm';
        
        numParameters = 1;
        
        inputText(numParameters).promptString = 'Select approach for second order IVA';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = char('newton', 'quasi');
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'second_order_opt_approach';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Enter max no. of iterations';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(1024);
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'maxIter';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Enter stopping tolerance';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = 1e-4;
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'termThreshold';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Enter learning rate';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = 0.2;
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'alpha0';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select stopping criteria';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = char('ChangeInCost', 'ChangeInW');
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'terminationCriterion';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Display statements';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = char('true', 'false');
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'verbose';
        inputText(numParameters).enable = 'on';
        
    case 'iva-l'
        %% IVA-L
        
        dlg_title = 'Select the options for the IVA-L algorithm';
        
        numParameters = 1;
        
        inputText(numParameters).promptString = 'Select weights initialization';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = char('gpca', 'random');
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'type';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Enter max no. of iterations';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(1024);
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'maxIter';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Enter stopping tolerance';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = 1e-4;
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'termThreshold';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Enter learning rate';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = 0.2;
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'alpha0';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select stopping criteria';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = char('ChangeInCost', 'ChangeInW');
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'terminationCriterion';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Display statements';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = char('true', 'false');
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'verbose';
        inputText(numParameters).enable = 'on';
        
        
    case 'complex infomax'
        %% Complex Infomax
        
        DEFAULT_LEARNING = 0.0001;
        
        % dialog Title
        dlg_title = 'Select the Options for the complex Infomax algorithm';
        
        numParameters = 1;
        
        inputText(numParameters).promptString = ' Select the nonlinearity parameter';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'1 - 2asin(u)', '1 - 2atan(u)', '1 - 2tanh(u)', '1 + 2acos(u)', 'Split tanh', 'Anem??ller' };
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'nonlinearity';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = ['Select Learning Level'];
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = num2str(DEFAULT_LEARNING);
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'learning';
        inputText(numParameters).enable = 'on';
        
        
    case 'complex fast ica'
        %% Complex Fast ICA
        
        
        % dialog Title
        dlg_title = 'Select the Options for the complex fastICA algorithm';
        
        numParameters = 1;
        
        
        inputText(numParameters).promptString = ' Select Symmetric or One by one estimation approach';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'Symmetric', 'One by one'};
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'complexestimation';
        inputText(numParameters).enable = 'on';
        
    case 'complex cmn'
        %% CMN algorithm
        
        % dialog Title
        dlg_title = 'Select the Options for the complex fastICA algorithm';
        
        numParameters = 1;
        
        
        inputText(numParameters).promptString = 'Select nonlinearity type';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'atanh','asinh', 'tanh', 'cosh', 'acosh', 'kurtosis', 'asin', 'tan' };
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'typestr';
        inputText(numParameters).enable = 'on';
        
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select whitening algorithm';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = { 'Normal', 'SUT'};
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'wtype';
        inputText(numParameters).enable = 'on';
        
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = 'Select data''s gaussian type';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = { 'subgaussian', 'supergaussian'};
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'subsuper';
        inputText(numParameters).enable = 'on';
        
        numParameters = numParameters + 1;
        
        inputText(numParameters).promptString = ' Select Symmetric or One by one estimation approach';
        inputText(numParameters).uiType = 'popup';
        inputText(numParameters).answerString = {'Symmetric', 'One by one'};
        inputText(numParameters).dataType = 'string';
        inputText(numParameters).tag = 'complexestimation';
        inputText(numParameters).enable = 'on';
        
end
%% End for ICA algorithm options

if exist('inputText', 'var')
    
    % Input dialog box
    answer = icatb_inputDialog('inputtext', inputText, 'Title', dlg_title, 'handle_visibility',  handle_visibility, 'windowStyle', windowStyle);
    
    
    % ICA options with flags and the values corresponding to it
    ICA_Options = cell(1, 2*length(answer));
    
    for i = 1:length(answer)
        ICA_Options{2*i - 1} = inputText(i).tag;
        ICA_Options{2*i} = answer{i};
    end
    
end