function icatb_enlwfc_lincopy(sesInfo)
    % Function supporting enl-WFC
    % Copies the explicitly nonlinear setup files 
    %    to a set of linear setup files
    % TReNDS 8/28/26 Cyrus Eierud
    
    %% Load defaults
    icatb_defaults;
    
    %% Enforce MAT file version
    global ENFORCE_MAT_FILE_VER;

    prefix_new = sesInfo.userInput.prefix;
    prefix_new = [prefix_new(1:length(prefix_new)-4) '-lin'];

    % Copy needed nii files
    copyfile([sesInfo.userInput.prefix 'Mask.nii'],[prefix_new 'Mask.nii']) 
    copyfile([sesInfo.userInput.prefix 'Mask_subsampled.nii'],[prefix_new 'Mask_subsampled.nii'])

    sesInfo.userInput.prefix = prefix_new;

    [pathstr, file_name, extn] = fileparts(sesInfo.userInput.maskFile);
    sesInfo.userInput.maskFile = [pathstr filesep prefix_new 'Mask' extn];
    
    clear pathstr file_name extn;
    [pathstr, file_name, extn] = fileparts(sesInfo.userInput.param_file);
    file_name_new = [prefix_new file_name(length(sesInfo.userInput.prefix)+1:end) extn];
    sesInfo.userInput.param_file = [pathstr filesep file_name_new];
    
    
    for n_sess = size(sesInfo.userInput.files,1)
       for n_sub = 1:size(sesInfo.userInput.files,2)
          clear pathstr file_name extn;
          [pathstr, file_name, extn] = fileparts(sesInfo.userInput.files(n_sess,n_sub).name);
          file_name_new = [prefix_new file_name(length(sesInfo.userInput.prefix)+1:end) extn];
          sesInfo.userInput.files(n_sess,n_sub).name = [pathstr filesep file_name_new];
       end
    end    
    
    sesInfo.userInput.lin_sidecar = false; % to stop coupying

    if (~isempty(ENFORCE_MAT_FILE_VER))
        save(sesInfo.userInput.param_file, 'sesInfo', ENFORCE_MAT_FILE_VER);
    else
        save(sesInfo.userInput.param_file, 'sesInfo');
    end        
    disp(['Saved file ', sesInfo.userInput.param_file, ' ...']);  

end

