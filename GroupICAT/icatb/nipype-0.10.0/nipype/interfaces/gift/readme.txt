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
matlab_cmd = '/data/mialab/users/srinivas/GIFT_Stand_alone/Linux_x86_64/GroupICATv4.0b_standalone_July_01_2020/run_groupica.sh /data/mialab/users/srinivas/GIFT_Stand_alone/Linux_x86_64/tmp_gica_stand_alone/v91/ '
gift.GICACommand.set_mlab_paths(matlab_cmd=matlab_cmd,use_mcr=True)
gc = gift.GICACommand()
gc.inputs.in_files = ['/data/mialab/users/srinivas/Example_Subjects/Visuomotor_data/sub01_vis/sub001.nii', '/data/mialab/users/srinivas/Example_Subjects/Visuomotor_data/sub02_vis/sub002.nii', '/data/mialab/users/srinivas/Example_Subjects/Visuomotor_data/sub03_vis/sub003.nii']
# first level pca components
gc.inputs.df = 30;
# Number of independent components
gc.inputs.dim = 16;
# Display results
gc.inputs.display_results = {'anatomical_file': '/data/mialab/users/srinivas/ch2bet.nii', 'slices_in_mm': list(range(-40, 76, 4)), 'convert_to_zscores':'yes', 'threshold':1.0, 'image_values':'positive', 'anatomical_plane':'axial'};
# Display network summary
gc.inputs.network_summary_opts = {'comp_network_names':{'VIS':[1,2,3], 'SM':[4,5],'AUD':6}, 'threshold' : 2.0, 'convert_to_z':'yes', 'format':'html'};
gc.run()

# Example 4: Evaluate any matlab command (only on the standard matlab path) or gift commands
from nipype.interfaces import gift   
matlab_cmd = '/data/mialab/users/srinivas/GIFT_Stand_alone/Linux_x86_64/GroupICATv4.0b_standalone_July30_2020/run_groupica.sh /data/mialab/users/srinivas/GIFT_Stand_alone/Linux_x86_64/tmp_gica_stand_alone/v91/ '
gift.evalGIFTCommand.set_mlab_paths(matlab_cmd=matlab_cmd,use_mcr=True)
ec = gift.evalGIFTCommand()
ec.inputs.file_name = '/data/mialab/users/srinivas/tmp_nipype/test_display.m';
ec.run();


# Example 5: Mancova example:

Stepwise regression:

from nipype.interfaces import gift   
matlab_cmd = '/data/mialab/users/srinivas/GIFT_Stand_alone/Linux_x86_64/GroupICATv4.0b_standalone_July30_2020/run_groupica.sh /data/mialab/users/srinivas/GIFT_Stand_alone/Linux_x86_64/tmp_gica_stand_alone/v91/ '
gift.MancovanCommand.set_mlab_paths(matlab_cmd=matlab_cmd,use_mcr=True)
mc = gift.MancovanCommand()
mc.inputs.ica_param_file = ['/path/to/ica_parameter_file']
mc.inputs.covariates = {'Age':['continuous', '/path/toage.txt', 'log'], 'Gender':['categorical', '/path/to/gender.txt']}
mc.inputs.comp_network_names = {'BG':21, 'VISUAL':[10, 12, 13]}
#mc.inputs.univariate_tests = {'Gender': ['Age'], 'Age': [], 'Gender_X_Age': []} #univariate tests examples.
#mc.inputs.univariate_tests={'Ttest':{'datasets':[ [i for i in range(1,26)] ], 'name':['HE']}} #One sample t-test
#mc.inputs.univariate_tests={'Ttest2':{'datasets':[[i for i in range(1,26)], [j for j in range(26, 51)]], 'name':['HE', 'SZ']}} #Two sample t-test
#mc.inputs.univariate_tests={'Ttest':{'datasets':[[i for i in range(1,26)], [j for j in range(26, 51)]], 'name':['Condition 1', 'Condition 2']}} #paired t-test
mc.inputs.numOfPCs = [4, 4, 4]
mc.inputs.TR = 2
mc.inputs.display = {'freq_limits':[0.1, 0.15], 'structFile':'/icatb_templates/ch2bet.nii', 't_threshold':1.0, 'image_values':'positive', 'threshdesc':'fdr', 'p_threshold':0.05};
mc.run()   

Univariate tests:

from nipype.interfaces import gift   
matlab_cmd = '/data/mialab/users/srinivas/GIFT_Stand_alone/Linux_x86_64/GroupICATv4.0b_standalone_July30_2020/run_groupica.sh /data/mialab/users/srinivas/GIFT_Stand_alone/Linux_x86_64/tmp_gica_stand_alone/v91/ '
gift.MancovanCommand.set_mlab_paths(matlab_cmd=matlab_cmd,use_mcr=True)
mc = gift.MancovanCommand()
mc.inputs.ica_param_file = ['/data/mialab/users/srinivas/Example_Subjects/mancova_sample_data/ica_output/rest_hcp_ica_parameter_info.mat']
mc.inputs.covariates = {'Age':['continuous', '/data/mialab/users/srinivas/Example_Subjects/mancova_sample_data/covariates/age.txt', 'log'], 'Gender':['categorical', '/data/mialab/users/srinivas/Example_Subjects/mancova_sample_data/covariates/gender.txt']}
mc.inputs.univariate_tests = {'Gender': ['Age'], 'Age': ['Gender'], 'Gender_X_Age': []} #univariate tests examples.
mc.inputs.comp_network_names = {'BG':21, 'VISUAL':[10, 12, 13]}
mc.inputs.TR = 2
mc.inputs.p_threshold = 0.2
mc.inputs.display = {'freq_limits':[0.1, 0.15], 't_threshold':1.0, 'image_values':'positive', 'threshdesc':'none', 'p_threshold':0.05};
mc.run()

One sample T-test example:

from nipype.interfaces import gift   
matlab_cmd = '/data/mialab/users/srinivas/GIFT_Stand_alone/Linux_x86_64/GroupICATv4.0b_standalone_July30_2020/run_groupica.sh /data/mialab/users/srinivas/GIFT_Stand_alone/Linux_x86_64/tmp_gica_stand_alone/v91/ '
gift.MancovanCommand.set_mlab_paths(matlab_cmd=matlab_cmd,use_mcr=True)
mc = gift.MancovanCommand()
mc.inputs.ica_param_file = ['/data/mialab/users/srinivas/Example_Subjects/mancova_sample_data/ica_output/rest_hcp_ica_parameter_info.mat']
mc.inputs.univariate_tests = {'Ttest': {'datasets':[ [i for i in range(1,51)] ], 'name':['HE']}} #data-sets should start with 1 (matlab syntax).
mc.inputs.comp_network_names = {'BG':21, 'VISUAL':[10, 12, 13]}
mc.inputs.TR = 2
mc.inputs.p_threshold = 0.2
mc.inputs.display = {'freq_limits':[0.1, 0.15], 't_threshold':1.0, 'image_values':'positive', 'threshdesc':'none', 'p_threshold':0.05};
mc.run()


Paired T-test example:

from nipype.interfaces import gift   
matlab_cmd = '/data/mialab/users/srinivas/GIFT_Stand_alone/Linux_x86_64/GroupICATv4.0b_standalone_July30_2020/run_groupica.sh /data/mialab/users/srinivas/GIFT_Stand_alone/Linux_x86_64/tmp_gica_stand_alone/v91/ '
gift.MancovanCommand.set_mlab_paths(matlab_cmd=matlab_cmd,use_mcr=True)
mc = gift.MancovanCommand()
mc.inputs.ica_param_file = ['/data/mialab/users/srinivas/Example_Subjects/mancova_sample_data/ica_output/rest_hcp_ica_parameter_info.mat']
mc.inputs.univariate_tests = {'Ttest': {'datasets':[ [i for i in range(1,26)], [i for i in range(26, 51)] ], 'name':['C1', 'C2']}} #data-sets should start with 1 (matlab syntax).
mc.inputs.comp_network_names = {'BG':21, 'VISUAL':[10, 12, 13]}
mc.inputs.TR = 2
mc.inputs.p_threshold = 0.2
mc.inputs.display = {'freq_limits':[0.1, 0.15], 't_threshold':1.0, 'image_values':'positive', 'threshdesc':'none', 'p_threshold':0.05};
mc.run()


Two sample T-test example:

from nipype.interfaces import gift   
matlab_cmd = '/data/mialab/users/srinivas/GIFT_Stand_alone/Linux_x86_64/GroupICATv4.0b_standalone_July30_2020/run_groupica.sh /data/mialab/users/srinivas/GIFT_Stand_alone/Linux_x86_64/tmp_gica_stand_alone/v91/ '
gift.MancovanCommand.set_mlab_paths(matlab_cmd=matlab_cmd,use_mcr=True)
mc = gift.MancovanCommand()
mc.inputs.ica_param_file = ['/data/mialab/users/srinivas/Example_Subjects/mancova_sample_data/ica_output/rest_hcp_ica_parameter_info.mat']
mc.inputs.univariate_tests = {'Ttest2': {'datasets':[ [i for i in range(1,26)], [i for i in range(26, 51)] ], 'name':['HE', 'SZ']}} #data-sets should start with 1 (matlab syntax).
mc.inputs.comp_network_names = {'BG':21, 'VISUAL':[10, 12, 13]}
mc.inputs.TR = 2
mc.inputs.p_threshold = 0.2
mc.inputs.display = {'freq_limits':[0.1, 0.15], 't_threshold':1.0, 'image_values':'positive', 'threshdesc':'none', 'p_threshold':0.05};
mc.run()

# Example 6: dfnc example

 import nipype.interfaces.gift as gift
 dc = gift.DFNCCommand()
 dc.inputs.ica_param_file = ['/path/to/ica_parameter_file']
 dc.inputs.comp_network_names = {'BG':21, 'VISUAL':[10, 12, 13]}
 dc.inputs.TR = 2
 dc.run()   

# Example 7: dfnc

from nipype.interfaces import gift   
matlab_cmd = '/data/mialab/users/srinivas/GIFT_Stand_alone/Linux_x86_64/GroupICATv4.0b_standalone_July30_2020/run_groupica.sh /data/mialab/users/srinivas/GIFT_Stand_alone/Linux_x86_64/tmp_gica_stand_alone/v91/ '
gift.DFNCCommand.set_mlab_paths(matlab_cmd=matlab_cmd,use_mcr=True)
dc = gift.DFNCCommand()
dc.inputs.ica_param_file = ['/data/mialab/users/srinivas/Example_Subjects/mancova_sample_data/ica_output/rest_hcp_ica_parameter_info.mat']
dc.inputs.comp_network_names = {'BG':21, 'VISUAL':[10, 12, 13]}
dc.inputs.TR = 2
dc.inputs.postprocess = {'num_clusters':5, 'display_results':1}
dc.run()
