import os
from nipype.interfaces.base import File, traits
from nipype.interfaces.matlab import MatlabInputSpec
from ..base import (isdefined, InputMultiPath)
from nipype.interfaces.gift.base import (GIFTCommand, GIFTCommandInputSpec)


class GICACommandInputSpec(GIFTCommandInputSpec): 
    """ Input specifications for groupica """
      
    in_files = InputMultiPath(File(exists=True), mandatory=True, desc="input file names (either single file name or a list)", sep=",")
    out_dir = traits.Str( mandatory = False, desc = 'Enter fullfile path of the results directory')
    mask = traits.Str( mandatory = False, desc = 'Enter file names using full path of the mask. If you wish to use default mask leave it empty')
    dim = traits.Int(mandatory = False, desc="dimensionality reduction into #num dimensions" "(default: automatic estimation)")	
    df = traits.Int(mandatory = False, desc="number of reduction steps used in the first pca step")	
    which_analysis = traits.Int(mandatory = False, desc="Options are 1, 2, and 3. 1 - standard group ica, 2 - ICASSO and 3 - MST. ")	
    group_ica_type = traits.Str(mandatory = False, desc="1 - Spatial ica, 2 - Temporal ica.")
    perfType = traits.Int(mandatory = False, desc="Options are 1, 2, and 3. 1 - maximize performance, 2 - less memory usage  and 3 - user specified settings. ")
    dummy_scans = traits.Int(mandatory = False, desc="enter dummy scans")
    prefix = traits.Str( mandatory = False, desc = 'Enter prefix to be appended with the output files')
    group_pca_type = traits.Str( mandatory = False, desc = 'options are subject specific and grand mean')
    backReconType = traits.Int( mandatory = False, desc = 'options are 1 - regular, 2 - spatial-temporal regression, 3 - gica3, 4 - gica, 5 - gig-ica')
    preproc_type = traits.Int( mandatory = False, desc = 'options are 1 - remove mean per timepoint, 2 - remove mean per voxel, 3 - intensity normalization, 4 - variance normalization')
    numReductionSteps = traits.Int( mandatory = False, desc = 'options are 1  and 2')
    doEstimation = traits.Int( mandatory = False, desc = 'options are 0 and 1')
    scaleType = traits.Int( mandatory = False, desc = 'options are 0 - No scaling, 1 - percent signal change, 2 - Z-scores')
    algoType = traits.Int( mandatory = False, desc = 'options are 1 - Infomax, 2 - Fast ica , ...')	
    refFiles = InputMultiPath(File(exists=True), argstr="-i %s", mandatory=False, position=0, desc="input file names (either single file name or a list)", sep=",")
    numWorkers = traits.Int( mandatory = False, desc = 'Number of parallel workers')	
    display_results = traits.Int( mandatory = False, desc = '0 - No display, 1 - HTML report, 2 - PDF')	
	

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
        	
        dummy_scans = 0;
        if isdefined(self.inputs.dummy_scans):
        		dummy_scans = self.inputs.dummy_scans;
        			
        batch_file_name = os.path.join(os.getcwd(), '%s_gica_batch.m' % (prefix)); 
        		
        commandstr = ["%% Batch script for running gica\n"];
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
        commandstr.append("%% Input mask\n");
        commandstr.append("maskFile = '%s';\n" % (mask));
        commandstr.append("\n");
        commandstr.append("dataSelectionMethod = 4;\n");
        commandstr.append("\n");
        commandstr.append("%% Input file patterns\n");
        commandstr.append("input_data_file_patterns = {");
        
        for n in range(len(self.inputs.in_files)):
        		commandstr.append("'%s';\n" % (self.inputs.in_files[n]));
        		  
        commandstr.append("};\n\n");
        commandstr.append("%% Dummy scans\n");
        commandstr.append("dummy_scans = %d;\n" % (dummy_scans));
        		
        if isdefined(self.inputs.group_pca_type):
        		commandstr.append("%% Group PCA type \n");
        		commandstr.append("group_pca_type = '%s';\n" % (self.inputs.group_pca_type));
        		
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
            display_results = self.inputs.display_results
        else:
            display_results = 0

        if isdefined(self.inputs.use_mcr) and self.inputs.use_mcr:
            display_results = 0
                
        commandstr.append("display_results = %d;\n" % (display_results));        
        
        fid = open(batch_file_name, "w+");
        fid.writelines(commandstr);
        fid.close();        
        script = "icatb_batch_file_run('%s')" % (batch_file_name)
        
        return script
		
class DFNCCommandInputSpec(GIFTCommandInputSpec):
    """ DFNC Inputs """
    
    ica_param_file = traits.Str( mandatory = True, desc = 'Enter fullfile path of the ICA parameter file')
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
    >>> dc.inputs.ica_param_file = /path/to/ica_parameter_file
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
        commandstr.append("ica_param_file = '%s';\n" % (self.inputs.ica_param_file));
	
        commandstr.append("%% TR in seconds \n");
        commandstr.append("TR = %f;\n" % (self.inputs.TR));
	
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
		
        if isdefined(self.inputs.use_mcr) and self.inputs.use_mcr:
            display_results = 0
	
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
    
    ica_param_file = traits.Str( mandatory = True, desc = 'Enter fullfile path of the ICA parameter file')
    out_dir = traits.Str( mandatory = False, desc = 'Enter fullfile path of the results directory')
    comp_network_names =  traits.Dict(mandatory=True, desc='dictionary containing network names and network values')
    TR =  traits.Float(mandatory=True, desc = "Enter experimental TR in seconds");
    features = traits.List(mandatory=False, desc = "Enter features like spatial maps, timecourses spectra, fnc correlations")
    covariates = traits.Dict(mandatory=True, desc='covariates. Each covariate must contain a list like category of covariate, file name, transformation (optional)')
    interactions = traits.List(mandatory=False, desc = "interaction terms")
    numOfPCs = traits.List(mandatory=True, desc = "Number of principal components for each feature")
    p_threshold =  traits.Float(mandatory=False, desc = "Enter p-threshold significance");
    feature_params = traits.Dict(mandatory=False, desc='Feature params')        
    
	
	
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
    >>> mc.inputs.ica_param_file = /path/to/ica_parameter_file
    >>> mc.inputs.comp_network_names = {'BG':21, 'VISUAL':[10, 12, 13]}
    >>> mc.inputs.numOfPCs = [4, 4, 4]
    >>> mc.inputs.covariates = {'Age':['continuous', '/path/toage.txt', 'log'], 'Gender':['categorical', '/path/to/gender.txt']}
    >>> mc.inputs.TR = 2
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
        commandstr.append("ica_param_file = '%s';\n" % (self.inputs.ica_param_file));
	
        commandstr.append("%% TR in seconds \n");
        commandstr.append("TR = %f;\n" % (self.inputs.TR));
	
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
        
        
        commandstr.append("%% Number of principal components for each feature \n");
        commandstr.append("numOfPCs = [");
        for n in range(len(self.inputs.numOfPCs)):
            commandstr.append("%d " % (self.inputs.numOfPCs[n]));     
        commandstr.append("];\n\n");    
            
            
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
                commandstr.append("%d " % (self.inputs.interactions[n]));
            commandstr.append("];\n\n"); 
        
        commandstr.append("%% Covariates \n");  
        commandstr.append("covariates = {");
        for keyN in self.inputs.covariates.keys():			
            commandstr.append("'%s', " % (keyN));			
            tmp_cov = self.inputs.covariates[keyN];		
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
            