How to Invoke GIFT on Nipype:

1. Download nipype-0.10.0 from http://nipype.readthedocs.io/en/latest/. Copy directory with sub-directories from GroupICATv4.0b/icatb/nipype-0.10.0 to nipype-0.10.0.
2. Use command pip install -e /path/to/local/nipype-0.10.0 at the command prompt.
3. Use ipython command to invoke python command prompt.

#Example 1: Single subject analysis

import nipype.interfaces.gift as gift
gc = gift.GICACommand()
# use full path
gc.inputs.in_files = '/path/to/swa.nii'
# Extract 20 components
gc.inputs.dim = 20
gc.run()  


#Example 2: Multi subject analysis

import nipype.interfaces.gift as gift
gc = gift.GICACommand()
gc.inputs.out_dir = '/path/to/output';
gc.inputs.numReductionSteps = 2
# use full path
gc.inputs.in_files = ['/path/to/swa1.nii', '/path/to/swa2.nii']
# Set 30 components in the first pca step
gc.inputs.df = 30
# Extract 20 components
gc.inputs.dim = 20
# set infomax algorithm
gc.inputs.algoType = 1
gc.run()  

#Example 3: Using MATLAB component runtime. Download standalone software from mialab.mrn.org/software/gift. 
# Optional: Set GICA_INSTALL_DIR environment variable using stand alone software path.

from nipype.interfaces import gift   
matlab_cmd = '/path/to/run_groupica.sh /path/to/compiler_runtime/v901/ '
gift.GICACommand.set_mlab_paths(matlab_cmd=matlab_cmd,use_mcr=True)
gc = gift.GICACommand()
gc.inputs.in_files = '/path/to/swa.nii'
# Use MDL estimation to detect the number of components.
gc.run()