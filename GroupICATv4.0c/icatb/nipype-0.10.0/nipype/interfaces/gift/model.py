import os
from nipype.interfaces.base import File, traits
from nipype.interfaces.matlab import MatlabInputSpec
from ..base import (isdefined, InputMultiPath)
from nipype.interfaces.gift.base import (GIFTCommand, GIFTCommandInputSpec)


class GICACommandInputSpec(GIFTCommandInputSpec): 
    """ Input specifications for groupica """
      
    in_files = InputMultiPath(File(exists=True), mandatory=True, desc="input file names (either single file name or a list)", sep=",")
    out_dir = traits.Str( mandatory = False, desc = 'Enter fullfile path of the results directory')
    mask = traits.Str( mandatory = False, desc = 'Options are default, average or enter fullfile name of the mask.')
    TR =  traits.List(mandatory=False, desc = "Enter experimental TR in seconds");
    dim = traits.Int(mandatory = False, desc="dimensionality reduction into #num dimensions" "(default: automatic estimation)")	
    df = traits.Int(mandatory = False, desc="number of reduction steps used in the first pca step")	
    which_analysis = traits.Int(mandatory = False, desc="Options are 1, 2, and 3. 1 - standard group ica, 2 - ICASSO and 3 - MST. ")	
    group_ica_type = traits.Str(mandatory = False, desc="1 - Spatial ica, 2 - Temporal ica.")
    perfType = traits.Int(mandatory = False, desc="Options are 1, 2, and 3. 1 - maximize performance, 2 - less memory usage  and 3 - user specified settings. ")
    dummy_scans = traits.List(mandatory = False, desc="enter dummy scans")
    prefix = traits.Str( mandatory = False, desc = 'Enter prefix to be appended with the output files')
    group_pca_type = traits.Str( mandatory = False, desc = 'options are subject specific and grand mean')
    pcaType = traits.Str(mandatory = False, desc = 'Options are Standard, Expectation Maximization, SVD, MPOWIT, STP'); 
    backReconType = traits.Int( mandatory = False, desc = 'options are 1 - regular, 2 - spatial-temporal regression, 3 - gica3, 4 - gica, 5 - gig-ica')
    preproc_type = traits.Int( mandatory = False, desc = 'options are 1 - remove mean per timepoint, 2 - remove mean per voxel, 3 - intensity normalization, 4 - variance normalization')
    numReductionSteps = traits.Int( mandatory = False, desc = 'options are 1  and 2')
    doEstimation = traits.Int( mandatory = False, desc = 'options are 0 and 1')
    scaleType = traits.Int( mandatory = False, desc = 'options are 0 - No scaling, 1 - percent signal change, 2 - Z-scores')
    algoType = traits.Int( mandatory = False, desc = 'options are 1 - Infomax, 2 - Fast ica , ...')	
    refFiles = InputMultiPath(File(exists=True), argstr="-i %s", mandatory=False, position=0, desc="input file names (either single file name or a list)", sep=",")
    numWorkers = traits.Int( mandatory = False, desc = 'Number of parallel workers')	
    #display_results = traits.Int( mandatory = False, desc = '0 - No display, 1 - HTML report, 2 - PDF')
    display_results = traits.Dict(mandatory=False, desc='dictionary containing results summary options')
    #network_summary_opts = {'comp_network_names':{'BG':17,'AUD':21,'SM':[23,24,29]}, 
    #'conn_threshold':0.1, 'structFile':'/path/to/ch2bet.nii', 'threshold':1};
    network_summary_opts = traits.Dict(mandatory=False, desc='dictionary containing network summary options')
    #icasso_opts = {'sel_mode':'randinit', 'num_ica_runs':10, 'min_cluster_size':8, 'max_cluster_size':10};
    icasso_opts = traits.Dict(mandatory=False, desc='dictionary containing icasso options')
    #mst_opts = {'num_ica_runs':10};
    mst_opts = traits.Dict(mandatory=False, desc='dictionary containing mst options')
    designMatrix = traits.List(mandatory = False, desc="enter SPM.mat file names in a list. You can enter one design for task based or equal to number of subjects for randomized designs")
    # example: ['right*bf(1)', 'left*bf(1)']
    regressors = traits.List(mandatory = False, desc="enter regressors from SPM.mat files")
	

class GICACommandOutputSpec(MatlabInputSpec):
    matlab_output = traits.Str( )

class GICACommand(GIFTCommand):
    """ Run Group ICA analysis using GIFT

    Returns
    -------

    matlab_output : capture of matlab output which may be
                    parsed by user to get computation results

    Examples
    --------

    >>> import nipype.interfaces.gift
    >>> gc = gift.GICACommand()
    >>> gc.inputs.in_files = '/path/to/swa.nii'
    >>> gc.run()    
    """
    input_spec = GICACommandInputSpec
    output_spec = GICACommandOutputSpec    

    def _make_matlab_command(self):
        """Implementation of GICA"""
        
        if isdefined(self.inputs.out_dir):
            os.chdir(self.inputs.out_dir);
            
        prefix = 'gica_cmd';
        if isdefined(self.inputs.prefix):
        		prefix = self.inputs.prefix;
                
        mask = '';
        if isdefined(self.inputs.mask):
        		mask = self.inputs.mask;    
        
        perfType = 1;
        if isdefined(self.inputs.perfType):
        		perfType = self.inputs.perfType;
                
        which_analysis = 1;
        if isdefined(self.inputs.which_analysis):
        		which_analysis = self.inputs.which_analysis;
        	
        icasso_opts = {'sel_mode':'randinit', 'num_ica_runs':10, 'min_cluster_size':8, 'max_cluster_size':10};
        
        mst_opts = {'num_ica_runs':10};
        
        if (isdefined(self.inputs.icasso_opts)):
            icasso_opts = self.inputs.icasso_opts;
            
        if (isdefined(self.inputs.mst_opts)):
            mst_opts = self.inputs.mst_opts;
        
        #Icasso options
        commandstr = ["%% Batch script for running gica\n"];
        commandstr.append("\n");
        commandstr.append("icasso_opts.sel_mode = '%s';\n" % (icasso_opts['sel_mode']));
        commandstr.append("icasso_opts.num_ica_runs = %d;\n" % (icasso_opts['num_ica_runs']));
        commandstr.append("icasso_opts.min_cluster_size = %d;\n" % (icasso_opts['min_cluster_size']));
        commandstr.append("icasso_opts.max_cluster_size = %d;\n" % (icasso_opts['max_cluster_size']));
        commandstr.append("\n");
        #MST options
        commandstr.append("mst_opts.num_ica_runs = %d;\n" % (mst_opts['num_ica_runs']));
        
        # Design info
        designMatrix = [];
        if isdefined(self.inputs.designMatrix):
        		designMatrix = self.inputs.designMatrix;
        
        keyword_design = 'no'
        if (len(designMatrix) > 0):
            if (len(designMatrix) == 1):
                keyword_design = 'same_sub_same_sess'
            else:
                keyword_design = 'diff_sub_diff_sess'
        
        commandstr.append("%% Design matrix/matrices\n");
        commandstr.append("keyword_designMatrix = '%s';\n" % (keyword_design));
        if (keyword_design == 'same_sub_same_sess'):
            commandstr.append("OnedesignMat = '%s';\n" % (designMatrix[0]));
        elif (keyword_design == 'diff_sub_diff_sess'):
             commandstr.append("input_design_matrices = {");
             if len(designMatrix) != len(self.inputs.in_files):
                 raise ValueError("Design matrix and input files must have the same list length")
             for nDesign in range(len(designMatrix)):
                 commandstr.append("'%s';\n" % (designMatrix[nDesign]))
             commandstr.append("};\n\n")
                 
        if isdefined(self.inputs.regressors):
            commandstr.append("%% Regressors\n");
            commandstr.append("refFunNames = {");
            for nRegress in range(len(self.inputs.regressors)):
                commandstr.append("'%s';\n" % (self.inputs.regressors[nRegress]))
            commandstr.append("};\n\n")
        
        dummy_scans = 0;
        if isdefined(self.inputs.dummy_scans):
        		dummy_scans = self.inputs.dummy_scans;
        			
        batch_file_name = os.path.join(os.getcwd(), '%s_gica_batch.m' % (prefix));
 		
        commandstr.append("\n");
        commandstr.append("%% Performance type\n");
        commandstr.append("perfType = %d;\n" % (perfType));
        commandstr.append("%% Reliability analysis\n");
        commandstr.append("which_analysis = %d;\n" % (which_analysis));
        commandstr.append("%% Output directory\n");
        commandstr.append("outputDir = '%s';\n" % (os.getcwd()));
        commandstr.append("\n");
        commandstr.append("%% Output files prefix\n");
        commandstr.append("prefix = '%s';\n" % (prefix));
        commandstr.append("\n");
        commandstr.append("dataSelectionMethod = 4;\n");
        commandstr.append("\n");
        commandstr.append("%% Input file patterns\n");
        commandstr.append("input_data_file_patterns = {");
        
        for n in range(len(self.inputs.in_files)):
        		commandstr.append("'%s';\n" % (self.inputs.in_files[n]));
        		  
        commandstr.append("};\n\n");
        commandstr.append("%% Dummy scans\n");
        commandstr.append("dummy_scans = %s;\n" % str(dummy_scans));
        
        if (mask.lower() == 'average'):
            commandstr.append("icatb_generateMask(input_data_file_patterns, 'outputdir', '%s', 'corr_threshold', 0.8);" % (os.getcwd()));
            mask = os.path.join(os.getcwd(), 'ica_analysisMask.nii');
            
        if (mask.lower() == 'default'):
            mask = '';
        
        commandstr.append("%% Input mask\n");
        commandstr.append("maskFile = '%s';\n" % (mask));
        commandstr.append("\n");
        
        if isdefined(self.inputs.TR):
            commandstr.append("%% TR in seconds \n");
            commandstr.append("TR = %s;\n" % str(self.inputs.TR));	
        
        if isdefined(self.inputs.group_pca_type):
        		commandstr.append("%% Group PCA type \n");
        		commandstr.append("group_pca_type = '%s';\n" % (self.inputs.group_pca_type));
                
                
        commandstr.append("%% PCA Algorithm\n");
        pcaType = 'Standard';
        if isdefined(self.inputs.pcaType):
        		pcaType = self.inputs.pcaType;       
        
        commandstr.append("pcaType = '%s';\n" % (pcaType));
        		
        backReconType = 4;
        if isdefined(self.inputs.backReconType):
        		backReconType = self.inputs.backReconType;
          
        commandstr.append("%% ICA Algorithm\n");
        algoType = 1;
        if isdefined(self.inputs.algoType):
        		algoType = self.inputs.algoType;   
          
        commandstr.append("algoType = %d;\n" % (algoType));          
        			
        commandstr.append("%% Back-reconstruction type\n");
        commandstr.append("backReconType = %d;\n" % (backReconType));
        		
        preproc_type = 4;
        	
        if isdefined(self.inputs.preproc_type):
        		preproc_type = self.inputs.preproc_type;
        
        commandstr.append("%% Pre-processing type\n");
        commandstr.append("preproc_type = %d;\n" % (preproc_type));
        			
        numReductionSteps = 2;
        
        if (len(self.inputs.in_files) > 2):
        		if isdefined(self.inputs.numReductionSteps):
        			numReductionSteps = self.inputs.numReductionSteps;
        		
        		
        if (numReductionSteps > 2):
        		numReductionSteps == 2;
        		
        commandstr.append("%% Number of data reduction steps\n");
        commandstr.append("numReductionSteps = %d;\n" % (numReductionSteps)); 
        	
        doEstimation = 0;			
        	
        if isdefined(self.inputs.doEstimation):	
        		doEstimation = self.inputs.doEstimation;
        	
        if not isdefined(self.inputs.dim):
        		doEstimation = 1;
        	
        commandstr.append("%% MDL Estimation \n");
        commandstr.append("doEstimation = %d;\n" % (doEstimation)); 
        	
        if (doEstimation == 0):
        		if (numReductionSteps == 1):
        			commandstr.append("%% Number of PC in the first PCA step\n");
        			commandstr.append("numOfPC1 = %d;\n" % (self.inputs.dim));
        				
        		if (numReductionSteps == 2):
        			dfValue = self.inputs.dim;
        			if isdefined(self.inputs.df):
        				dfValue = self.inputs.df;
        			commandstr.append("%% Number of PC in the first PCA step\n");
        			commandstr.append("numOfPC1 = %d;\n" % (dfValue));
        			commandstr.append("%% Number of PC in the second PCA step\n");
        			commandstr.append("numOfPC2 = %d;\n" % (self.inputs.dim));
        
        scaleType = 2;
        if isdefined(self.inputs.scaleType):
        		scaleType = self.inputs.scaleType;		
        commandstr.append("%% Scaling type \n");
        commandstr.append("scaleType = %d;\n" % (scaleType));		
        	
        	
        if isdefined(self.inputs.refFiles):	
        		commandstr.append("%% Spatial references \n");
        		commandstr.append("refFiles = {");
        		for n in range(len(self.inputs.refFiles)):			
        			commandstr.append("'%s';\n" % (self.inputs.refFiles[n]));
        		commandstr.append("};\n\n");
        		
        if isdefined(self.inputs.numWorkers):
        		commandstr.append("%% Parallel info \n");
        		commandstr.append("parallel_info.mode='parallel';");
        		commandstr.append("parallel_info.num_workers = %d;" % (self.inputs.numWorkers));

        commandstr.append("%% Report generator \n");                
        if isdefined(self.inputs.display_results):            
            tmp_disp_results = self.inputs.display_results;
        else:
            tmp_disp_results = 0;
        
        display_results = tmp_disp_results;        
        if not ((type(display_results) == int) or  (type(display_results) == str)):
            for keyD in display_results.keys():
                tmp_disp_results = display_results[keyD];
                if (type(tmp_disp_results) == str):
                    commandstr.append("display_results.{} = '{}';\n".format(keyD, tmp_disp_results));
                else:
                    commandstr.append("display_results.{} = {};\n".format(keyD, tmp_disp_results));
        else:
            if (type(display_results) == str):
                commandstr.append("display_results = '{}';\n".format(display_results));
            else:
                commandstr.append("display_results = {};\n".format(display_results));
        
        # network summary options
        if isdefined(self.inputs.network_summary_opts):
			network_summary_opts = self.inputs.network_summary_opts;
            		commandstr.append("%% Network summary options \n");
            		comp_network_names = network_summary_opts['comp_network_names'];
            		commandstr.append("network_summary_opts.comp_network_names = {");
            		for keyN in comp_network_names.keys():
                		commandstr.append("'%s', " % (keyN));
                		if type(comp_network_names[keyN]) == int:
                    			vals = [comp_network_names[keyN]];
                		else:
                    			vals = comp_network_names[keyN];
                		commandstr.append("%s;\n" % str(vals));
                	commandstr.append("};\n\n");

                	# Network summary threshold
                	try:
                    		network_threshold = self.inputs.network_summary_opts['threshold'];
                	except:
                    		network_threshold = 1;
                	commandstr.append("network_summary_opts.threshold = %s;\n" % str(network_threshold));

                    # Connectivity threshold
                	try:
                            conn_threshold = self.inputs.network_summary_opts['conn_threshold'];
                            commandstr.append("network_summary_opts.conn_threshold = %s;\n" % str(conn_threshold));
			except: 
			    pass;
                    # Anatomical file
                	try:
                    		structFile = self.inputs.network_summary_opts['structFile'];
                    		commandstr.append("network_summary_opts.structFile = '%s';\n" % (structFile));
                	except:
                            pass;
                    
                	commandstr.append("network_summary_opts.save_info = 1;\n");

                	# network summary format
                	try:
                    		file_format = self.inputs.network_summary_opts['format'];
                	except:
                    		file_format = "html";

                	commandstr.append("network_summary_opts.format = '%s';\n" % (file_format));

                	# Convert to z-scores
                	try:
                    		convert_to_z = self.inputs.network_summary_opts['convert_to_z'];
                	except:
                    		convert_to_z = "yes";

	                commandstr.append("network_summary_opts.convert_to_z = '%s';\n" % (convert_to_z));
        
        fid = open(batch_file_name, "w+");
        fid.writelines(commandstr);
        fid.close();        
        script = "icatb_batch_file_run('%s')" % (batch_file_name)
        
        return script    


class evalGIFTCommandInputSpec(GIFTCommandInputSpec):
    """ Input specifications for the script """

    file_name = traits.Str( mandatory = True, desc = 'Enter full file name.');

class evalGIFTCommandOutputSpec( MatlabInputSpec):
    
    matlab_output = traits.Str( )

class evalGIFTCommand(GIFTCommand):
    """ Evaluate any matlab command """

    input_spec = evalGIFTCommandInputSpec
    output_spec = evalGIFTCommandOutputSpec

    def _make_matlab_command(self):
        """Implementation of eval script"""
        
        script = "icatb_eval_script('%s')" % (self.inputs.file_name)
        
        return script;
    
    
		
class DFNCCommandInputSpec(GIFTCommandInputSpec):
    """ DFNC Inputs """
    
    ica_param_file = traits.List( mandatory = True, desc = 'Enter fullfile path of the ICA parameter file in a list')
    out_dir = traits.Str( mandatory = False, desc = 'Enter fullfile path of the results directory')
    comp_network_names =  traits.Dict(mandatory=True, desc='dictionary containing network names and network values')
    TR =  traits.Float(mandatory=True, desc = "Enter experimental TR in seconds");
    dfnc_params =  traits.Dict(mandatory=False, desc='dictionary containing dfnc parameters')
    postprocess = traits.Dict(mandatory=False, desc='dictionary containing post-processing parameters')
	
	
class DFNCCommandOutputSpec( MatlabInputSpec):
	matlab_output = traits.Str( )	
	
	
class DFNCCommand(GIFTCommand):
    """ Run dFNC using GIFT

    Returns
    -------

    matlab_output : capture of matlab output which may be
                    parsed by user to get computation results

    Examples
    --------

    >>> import nipype.interfaces.gift
    >>> dc = gift.DFNCCommand()
    >>> dc.inputs.ica_param_file = ['/path/to/ica_parameter_file']
    >>> dc.inputs.comp_network_names = {'BG':21, 'VISUAL':[10, 12, 13]}
    >>> dc.inputs.TR = 2
    >>> dc.run()   
    """
    input_spec = DFNCCommandInputSpec
    output_spec = DFNCCommandOutputSpec    

    def _make_matlab_command(self):
        """Implementation of dFNC """
        
        if isdefined(self.inputs.out_dir):
            os.chdir(self.inputs.out_dir);			
            
        batch_file_name = os.path.join(os.getcwd(), 'nipype_dfnc_batch.m'); 
	
        commandstr = ["%% Batch script for running dFNC\n"];
        commandstr.append("\n");
        commandstr.append("outputDir = '%s';\n" % (os.getcwd()));
	
        commandstr.append("%% ICA parameter file name \n");
        commandstr.append("ica_param_file = {");
        for nlist in self.inputs.ica_param_file:
            commandstr.append("'%s'  " % (nlist));
        commandstr.append("};\n");
        
        commandstr.append("%% TR in seconds \n");
        commandstr.append("TR = %s;\n" % str(self.inputs.TR));
	
        commandstr.append("%% Network names and values \n");	
        commandstr.append("comp_network_names = {");
	
        for keyN in self.inputs.comp_network_names.keys():			
            commandstr.append("'%s', [" % (keyN));			
            if type(self.inputs.comp_network_names[keyN]) == int:
                vals = [self.inputs.comp_network_names[keyN]];
            else:
                vals = self.inputs.comp_network_names[keyN];			
                
            for n in range(len(vals)):
                commandstr.append("%d " % (vals[n]));
            commandstr.append("];\n");
        commandstr.append("};\n\n");
			
        try:
            tc_detrend = self.inputs.dfnc_params['tc_detrend'];		
        except:
            tc_detrend = 3;
	
        try:
            tc_despike = self.inputs.dfnc_params['tc_despike'];
        except:
            tc_despike = 'yes';
		
        try:
            tc_filter = self.inputs.dfnc_params['tc_filter'];
        except:
            tc_filter = 0.15;
			
        try:
            method = self.inputs.dfnc_params['method'];
        except:
            method = 'none';
		
	
        try:
            wsize = self.inputs.dfnc_params['wsize'];
        except:
            wsize = 30;
	
        try:
            window_alpha = self.inputs.dfnc_params['window_alpha'];
        except:
           window_alpha = 3;
		
	
        try:
            num_repetitions = self.inputs.dfnc_params['num_repetitions'];	
        except:
            num_repetitions = 10;
		
		
        commandstr.append("%% dfnc parameters \n");
        commandstr.append("dfnc_params.tc_detrend = %d; \n" % (tc_detrend));
        commandstr.append("dfnc_params.tc_despike = '%s'; \n" % (tc_despike));
        commandstr.append("dfnc_params.tc_despike = %d; \n" % (tc_filter));
        commandstr.append("dfnc_params.method = '%s'; \n" % (method));
        commandstr.append("dfnc_params.wsize = %d; \n" % (wsize));
        commandstr.append("dfnc_params.window_alpha = %d; \n" % (window_alpha));
        commandstr.append("dfnc_params.num_repetitions = %d; \n" % (num_repetitions));
	
        try:
            self.inputs.dfnc_params['filesList']
            commandstr.append("dfnc_params.tc_covariates.filesList = {");
            for n in range(len(self.inputs.dfnc_params['filesList'])):			
                commandstr.append("'%s';\n" % (self.inputs.dfnc_params['filesList'][n]));
            commandstr.append("};\n\n");
        except:
            commandstr.append("\n");
	
        try:
            num_clusters = self.inputs.postprocess['num_clusters']
        except:
            num_clusters = 3
		
        try:
            ica_comps = self.inputs.postprocess['ica_comps']
        except:
            ica_comps = 3
			
        try:
            ica_algorithm = self.inputs.postprocess['ica_algorithm']
        except:
            ica_algorithm = 1
		
        try:
            num_ica_runs = self.inputs.postprocess['num_ica_runs']
        except:
            num_ica_runs = 5
			
        try:
            dmethod = self.inputs.postprocess['dmethod']
        except:
            dmethod = 'city'
			
        try:
            kmeans_max_iter = self.inputs.postprocess['kmeans_max_iter']
        except:
            kmeans_max_iter = 150
			
        try:
            display_results = self.inputs.postprocess['display_results']
        except:
            display_results = 1
		
        #if isdefined(self.inputs.use_mcr) and self.inputs.use_mcr:
        #    display_results = 0
	
        try:
            regressCovFile = self.inputs.postprocess['regressCovFile'];
        except:
            regressCovFile = '';

        commandstr.append("%% Postprocess parameters \n");
        commandstr.append("postprocess.num_clusters = %d; \n" % (num_clusters));
        commandstr.append("postprocess.ica.num_comps = %d; \n" % (ica_comps));
        commandstr.append("postprocess.ica.algorithm = %d; \n" % (ica_algorithm));
        commandstr.append("postprocess.regressCovFile = '%s'; \n" % (regressCovFile));
        commandstr.append("postprocess.ica.num_ica_runs = '%s'; \n" % (num_ica_runs));
        commandstr.append("postprocess.dmethod = '%s'; \n" % (dmethod));
        commandstr.append("postprocess.kmeans_max_iter = %d; \n" % (kmeans_max_iter));
        commandstr.append("postprocess.display_results = %d; \n" % (display_results));	
        
        fid = open(batch_file_name, "w+");
        fid.writelines(commandstr);
        fid.close();
        script = "icatb_dfnc_batch('%s')" % (batch_file_name)
        return script


class MancovanCommandInputSpec(GIFTCommandInputSpec):
    """ Mancovan Inputs """
    
    ica_param_file = traits.List( mandatory = True, desc = 'Enter fullfile path of the ICA parameter file in a list')
    out_dir = traits.Str( mandatory = False, desc = 'Enter fullfile path of the results directory')
    comp_network_names =  traits.Dict(mandatory=False, desc='dictionary containing network names and network values')
    TR =  traits.Float(mandatory=False, desc = "Enter experimental TR in seconds");
    features = traits.List(mandatory=False, desc = "Enter features like spatial maps, timecourses spectra, fnc correlations")
    covariates = traits.Dict(mandatory=False, desc='covariates. Each covariate must contain a list like category of covariate, file name, transformation (optional)')
    interactions = traits.List(mandatory=False, desc = "interaction terms")
    numOfPCs = traits.List(mandatory=False, desc = "Number of principal components for each feature")
    p_threshold =  traits.Float(mandatory=False, desc = "Enter p-threshold significance");
    feature_params = traits.Dict(mandatory=False, desc='Feature params');
    write_stats_info = traits.Int( mandatory = False, desc = '1 - Write stats info which is used later in aggregating results across sites.')
    univariate_tests = traits.Dict(mandatory=False, desc='Univariate tests to specify');
    display = traits.Dict(mandatory=False, desc='Display');        
    
	
	
class MancovanCommandOutputSpec( MatlabInputSpec):
	matlab_output = traits.Str( )	
	
	
class MancovanCommand(GIFTCommand):
    """ Run Mancovan using GIFT

    Returns
    -------

    matlab_output : capture of matlab output which may be
                    parsed by user to get computation results

    Examples
    --------

    >>> import nipype.interfaces.gift
    >>> mc = gift.MancovanCommand()
    >>> mc.inputs.ica_param_file = ['/path/to/ica_parameter_file'];
    >>> mc.inputs.feature_params = {'spatial maps': {'sm_center': 'No', 'sm_mask': ''},
                          'spectra': {'spectra_detrend' : 3, 'spectra_normalize_subs' : 'yes', 'spectra_transform' : 'yes'},
                          'fnc': {'fnc_tc_detrend' : 3, 'fnc_tc_despike' : 'yes', 'fnc_tc_filter' : 0.15}};
    >>> #mc.inputs.univariate_tests = {'Gender': ['Age'], 'Age': [], 'Gender_X_Age': []} #univariate tests examples.
    >>> #mc.inputs.univariate_tests={'Ttest':{'datasets':[ [i for i in range(1,26)] ], 'name':['HE']}} #One sample t-test
    >>> #mc.inputs.univariate_tests={'Ttest2':{'datasets':[[i for i in range(1,26)], [j for j in range(26, 51)]], 'name':['HE', 'SZ']}} #Two sample t-test
    >>> #mc.inputs.univariate_tests={'Ttest':{'datasets':[[i for i in range(1,26)], [j for j in range(26, 51)]], 'name':['Condition 1', 'Condition 2']}} #paired t-test
    >>> mc.inputs.comp_network_names = {'BG':21, 'VISUAL':[10, 12, 13]}
    >>> mc.inputs.numOfPCs = [4, 4, 4]
    >>> mc.inputs.covariates = {'Age':['continuous', '/path/toage.txt', 'log'], 'Gender':['categorical', '/path/to/gender.txt']}
    >>> mc.inputs.TR = 2
    >>> mc.inputs.feature_params = {'spatial maps': {'sm_center': 'No', 'sm_mask': '/tmp/myMask.nii'},
                          'spectra': {'spectra_detrend' : 3, 'spectra_normalize_subs' : 'yes', 'spectra_transform' : 'yes'},
                          'fnc': {'fnc_tc_detrend' : 3, 'fnc_tc_despike' : 'yes', 'fnc_tc_filter' : 0.15}};
    >>> mc.inputs.display = {'freq_limits':[0.1, 0.15], 'structFile':'/icatb_templates/ch2bet.nii', 
    't_threshold':1.0, 'image_values':'positive', 'threshdesc':'fdr', 'p_threshold':0.05, 'display_connectogram': 1, 
    'compFiles':'/icatb_templates/neuromark_53.nii'};
    >>> mc.run()   
    """
    input_spec = MancovanCommandInputSpec
    output_spec = MancovanCommandOutputSpec  
    
    def _make_matlab_command(self):
        """Implementation of Mancovan """
        
        if isdefined(self.inputs.out_dir):
            os.chdir(self.inputs.out_dir);			
            
        batch_file_name = os.path.join(os.getcwd(), 'nipype_mancovan_batch.m'); 
        
        commandstr = ["%% Batch script for running Mancovan\n"];
        commandstr.append("\n");
        commandstr.append("outputDir = '%s';\n" % (os.getcwd()));
	
        commandstr.append("%% ICA parameter file name \n");
        commandstr.append("ica_param_file = {");
        for nlist in self.inputs.ica_param_file:
            commandstr.append("'%s'  " % (nlist));
        commandstr.append("};\n");
        
        if (isdefined(self.inputs.display)):
            commandstr.append("%% display parameters \n");
            for keyN in self.inputs.display.keys():
                commandstr.append("display.{} = [".format(keyN));
                if type(self.inputs.display[keyN]) == list:
                    vals = self.inputs.display[keyN];
                elif type(self.inputs.display[keyN]) == str:
                    vals = ["'{}'".format(self.inputs.display[keyN])];
                else:
                    vals = [self.inputs.display[keyN]];
                    
                for n in range(len(vals)):
                    commandstr.append("{} ".format(vals[n]));
                commandstr.append("];");
                commandstr.append("\n");        
    
        # return the params defined for aggregating mancova stats
        if (self.inputs.ica_param_file[0].endswith('stats_info.mat')):
            fid = open(batch_file_name, "w+");
            fid.writelines(commandstr);
            fid.close();
            script = "icatb_mancovan_batch('%s')" % (batch_file_name)
            return script
        
        # Feature parameters
        feature_params = {'spatial maps': {'sm_center': 'No', 'sm_mask': ''},
                          'spectra': {'spectra_detrend' : 3, 'spectra_normalize_subs' : 'yes', 'spectra_transform' : 'yes'},
                          'fnc': {'fnc_tc_detrend' : 3, 'fnc_tc_despike' : 'yes', 'fnc_tc_filter' : 0.15}};
        
        if isdefined(self.inputs.feature_params):
            for nFeat in self.inputs.feature_params:
                current_feature = self.inputs.feature_params[nFeat];
                for nCFeat in current_feature:
                    feature_params[nFeat][nCFeat] = current_feature[nCFeat];
        
        for nFeat in feature_params.keys():
            current_feature = feature_params[nFeat];
            for nCFeat in current_feature.keys():
                tmp = feature_params[nFeat][nCFeat];
                if type(tmp) == int:
                    commandstr.append("feature_params.%s = %d;\n" % (nCFeat, tmp));
                elif type(tmp) == float:
                    commandstr.append("feature_params.%s = %f;\n" % (nCFeat, tmp));
                else:
                    commandstr.append("feature_params.%s = '%s';\n" % (nCFeat, tmp));
        
        
        write_stats_info = 1
        if (isdefined(self.inputs.write_stats_info) and self.inputs.write_stats_info is not None):
            write_stats_info = self.inputs.write_stats_info;
        
        commandstr.append("\n")
        commandstr.append("%% Write stats info \n");
        commandstr.append("write_stats_info = %s;\n" % str(write_stats_info));
        
        commandstr.append("%% TR in seconds \n");
        commandstr.append("TR = %s;\n" % str(self.inputs.TR));
	
        commandstr.append("%% Network names and values \n");	
        commandstr.append("comp_network_names = {");
	
        for keyN in self.inputs.comp_network_names.keys():			
            commandstr.append("'%s', [" % (keyN));			
            if type(self.inputs.comp_network_names[keyN]) == int:
                vals = [self.inputs.comp_network_names[keyN]];
            else:
                vals = self.inputs.comp_network_names[keyN];			
                
            for n in range(len(vals)):
                commandstr.append("%d " % (vals[n]));
            commandstr.append("];\n");
        commandstr.append("};\n\n");
            

        commandstr.append("%% Features \n");            
        if  isdefined(self.inputs.features):
            features = self.inputs.features;
        else:
            features = ['spatial maps', 'timecourses spectra', 'fnc correlations'];
			
        commandstr.append("features = {");
        for n in range(len(features)):
            commandstr.append("'%s';" % (features[n]));     
        commandstr.append("};\n\n");
        
        
        if  isdefined(self.inputs.numOfPCs):
            commandstr.append("%% Number of principal components for each feature \n");
            commandstr.append("numOfPCs = [");
            for n in range(len(self.inputs.numOfPCs)):
                commandstr.append("%d " % (self.inputs.numOfPCs[n]));
            commandstr.append("];\n\n");    
        
        
        if isdefined (self.inputs.univariate_tests):
            univariate_tests = self.inputs.univariate_tests;
            univ_keys = list(univariate_tests.keys());
            if (univ_keys[0].lower() == 'ttest2' or univ_keys[0].lower() == 'ttest'):
                
                ttest_datasets = univariate_tests[univ_keys[0]]['datasets'];
                
                if (len(ttest_datasets) > 2):
                    raise Exception('Dataset error: Nested lists cannot be greater than length of 2');
                
                try:
                    ttest_group_names = univariate_tests[univ_keys[0]]['name'];
                except:
                    ttest_group_names=['Group'];
                    if (univ_keys[0].lower() == 'ttest2'):
                        ttest_group_names = ['Group 1', 'Group 2'];
                    else:
                        if (len(ttest_datasets) == 2):
                            ttest_group_names = ['Condition 1', 'Condition 2'];
                
                commandstr.append("%% Univariate tests \n");
                commandstr.append("univariate_tests = {'%s', {" % univ_keys[0]);
                for nD in ttest_datasets:
                    commandstr.append("%s " % str(nD));
                
                commandstr.append("}, {");
                for nD in ttest_group_names:
                    commandstr.append("'%s' " % nD);
                
                commandstr.append("}};\n\n");
                
            else:
                commandstr.append("%% Univariate tests \n");
                commandstr.append("univariate_tests = {");
                for nD in univariate_tests.keys():
                    commandstr.append("'%s', {" % nD);
                    for nStr in univariate_tests[nD]:
                        commandstr.append("'%s' " % nStr);
                    commandstr.append("};\n");
                commandstr.append("};\n\n");
        
        if isdefined(self.inputs.p_threshold):
            p_threshold = self.inputs.p_threshold
        else:
            p_threshold = 0.01;            
        
        commandstr.append("%% p-threshold \n");    
        commandstr.append("p_threshold = %f;\n" % (p_threshold));
                    
        if isdefined(self.inputs.interactions):
            commandstr.append("%% Interaction terms if any \n");
            commandstr.append("interactions = [");
            for n in range(len(self.inputs.interactions)):
                commandstr.append("%s; " % str(self.inputs.interactions[n]));
            commandstr.append("];\n\n"); 
        
        if not isdefined(self.inputs.covariates):
            covariates = {};
        else:
            covariates = self.inputs.covariates;
        commandstr.append("%% Covariates \n");  
        commandstr.append("covariates = {");
        for keyN in covariates.keys():			
            commandstr.append("'%s', " % (keyN));			
            tmp_cov = covariates[keyN];		
            commandstr.append("'%s', " % (tmp_cov[0]));	
            commandstr.append("'%s', " % (tmp_cov[1]));  
            try:
                tmp_transform = tmp_cov[2];
            except:
                tmp_transform = '';
            commandstr.append("'%s';\n" % (tmp_transform)); 
        commandstr.append("};\n\n");               
                    
        fid = open(batch_file_name, "w+");
        fid.writelines(commandstr);
        fid.close();
        script = "icatb_mancovan_batch('%s')" % (batch_file_name)
        return script            
            
