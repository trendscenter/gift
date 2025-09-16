function icatb_compile_mex_files
%% Function to compile mex binaries of GroupICAT
%
%

compileGlassoFile;

compileSPMFiles;

oldDir = pwd;

%% Group ICA MEX directory
gift_mex_dir = fullfile(fileparts(which('groupica.m')), 'icatb_mex_files', 'src');

cd(gift_mex_dir);

%% Get mex extension
mexExtension = mexext;

%% Matlab root directory
matlabRootDir = matlabroot;

%% Matlab version
matlab_version = icatb_get_matlab_version;

MWSIZE_MACRO = '';
largeArrayDims = '';
if (matlab_version >= 2008)
    MWSIZE_MACRO = '-DMWSIZE';
end

%% Operating system bit
OSBIT = 32;
if (strcmp(mexExtension(end-1:end), '64'))
    OSBIT = 64;
end

%% SINGLE PRECISION MACRO
SINGLE_SUPPORT_MACRO = '';
if (matlab_version > 13)
    if (OSBIT == 32)
        %% 32 bit
        SINGLE_SUPPORT_MACRO = '-DSINGLE_SUPPORT';
    else
        %% 64 bit
        if (matlab_version >= 2008)
            SINGLE_SUPPORT_MACRO = '-DSINGLE_SUPPORT';
            largeArrayDims = '-largeArrayDims';
        end
    end
end

%% Mex files
files = {'icatb_eig_symm_all.c', 'icatb_eig_symm_sel.c'};

%% PC MACRO
PC_MACRO = '';

if (ispc)
    %% Windows
    
    PC_MACRO = '-DPC';
    OSDir = ['win', num2str(OSBIT)];
    
    %% Lapack and Compiler options files
    if (OSBIT == 32)
        lapack_lib = fullfile(matlabRootDir, 'extern', 'lib', OSDir, 'lcc', 'libmwlapack.lib');
        compiler_bat = fullfile(matlabRootDir, 'bin', OSDir, 'mexopts', 'lccopts.bat');
    else
        lapack_lib = fullfile(matlabRootDir, 'extern', 'lib', OSDir, 'microsoft', 'libmwlapack.lib');
        compiler_bat = '';
    end
    
end

%% Loop over files
for nF = 1:length(files)
    fprintf('\n');
    current_file = files{nF};
    [p, fName] = fileparts(current_file);
    disp(['Compiling ', current_file, ' ...']);
    if (ispc)
        if (OSBIT == 32)
            commandToRun = ['mex "', current_file, '" ', largeArrayDims, ' ', PC_MACRO, ' ', SINGLE_SUPPORT_MACRO, ' ', MWSIZE_MACRO, ' "', lapack_lib, '" -f "', compiler_bat, '"'];
        else
            commandToRun = ['mex "', current_file, '" ', largeArrayDims, ' ', PC_MACRO, ' ', SINGLE_SUPPORT_MACRO, ' ', MWSIZE_MACRO, ' "', lapack_lib, '"'];
        end
    else
        commandToRun = ['mex "', current_file, '" ', largeArrayDims, ' ', PC_MACRO, ' ', SINGLE_SUPPORT_MACRO, ' ', MWSIZE_MACRO, ' -lmwlapack'];
    end
    disp(['Executing ', commandToRun]);
    eval(commandToRun);
    movefile([fName, '.', mexExtension], '..');
end
%% End of loop over files
fprintf('\nDone\n');


cd(oldDir);

function compileSPMFiles
%% Compile SPM files
%

if (~ispc)
    sourceDir = fullfile(fileparts(which('groupica.m')), 'icatb_spm8_files', 'src');
    commandToRun = ['cd "', sourceDir, '";make all install'];
    [s, m] = system(commandToRun);
    if (s)
        fprintf('\n');
        disp(m);
        fprintf('\n Could not compile SPM mex files. Try copying MEX files manually from spm8 to icatb/icatb_spm8_files. See group ICA manual \n\n');
    end
end

function compileGlassoFile
%% Compile Glasso file
%
outDir = pwd;
%if (~ispc)
try
    sourceDir = fullfile(fileparts(which('groupica.m')), 'toolbox', 'Graphical_Lasso');
    cd(sourceDir);
    mex glasso.F;
catch
    disp(['!!!Problem creating ', fullfile(sourceDir, ['glasso.', mexext]), ' file. Manually compile ', fullfile(sourceDir, 'glasso.F')]);
end
%end
cd(outDir);