# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""The GIFT module provides basic functions for interfacing with some of the GIFT  tools.

In order to use the standalone MCR version of GIFT, you need to ensure that
the following commands are executed at the beginning of your script::

   from nipype.interfaces import gift   
   matlab_cmd = '/path/to/run_groupica.sh /path/to/compiler_runtime/v901/ '
   gift.GICACommand.set_mlab_paths(matlab_cmd=matlab_cmd,use_mcr=True)
"""

__docformat__ = 'restructuredtext'

# Standard library imports
import os

# Local imports
from ..base import (BaseInterface, traits, isdefined, InputMultiPath,
                    BaseInterfaceInputSpec, Directory, Undefined)
from ..matlab import MatlabCommand

class GIFTCommandInputSpec(BaseInterfaceInputSpec):
    matlab_cmd = traits.Str(desc='matlab command to use')
    paths = InputMultiPath(Directory(), desc='Paths to add to matlabpath')
    mfile = traits.Bool(True, desc='Run m-code using m-file', usedefault=True)
    use_mcr = traits.Bool(desc='Run m-code using GIFT MCR')     
	
class GIFTCommandOutputSpec( BaseInterfaceInputSpec):
    matlab_output = traits.Str( )	

class GIFTCommand(BaseInterface):
    """Extends `BaseInterface` class to implement GIFT specific interfaces.

    WARNING: Pseudo prototype class, meant to be subclassed
    """
    input_spec = GIFTCommandInputSpec
    output_spec = GIFTCommandOutputSpec
    
    _matlab_cmd = None
    _paths = None
    _use_mcr = None

    def __init__(self, **inputs):
        super(GIFTCommand, self).__init__(**inputs)
        self.inputs.on_trait_change(self._matlab_cmd_update, ['matlab_cmd','mfile','paths','use_mcr'])
        self._find_mlab_cmd_defaults()
        self._check_mlab_inputs()
        self._matlab_cmd_update()

    @classmethod
    def set_mlab_paths(cls, matlab_cmd=None, paths=None, use_mcr=None):
        cls._matlab_cmd = matlab_cmd
        cls._paths = paths
        cls._use_mcr = use_mcr

    def _find_mlab_cmd_defaults(self):
        # check if the user has set environment variables to enforce
        # the standalone (MCR) version of GIFT        
        if self._use_mcr:
            self._use_mcr = True
          

    def _matlab_cmd_update(self):
        # MatlabCommand has to be created here,
        # because matlab_cmb is not a proper input
        # and can be set only during init	
        matlab_cmd_str = self.inputs.matlab_cmd	
        if isdefined(self.inputs.use_mcr) and self.inputs.use_mcr:
            if not matlab_cmd_str[-1] == " ":
                matlab_cmd_str = matlab_cmd_str + " "
        self.mlab = MatlabCommand(matlab_cmd=matlab_cmd_str,
                                  mfile=self.inputs.mfile,
                                  paths=self.inputs.paths)       
        self.mlab.inputs.script_file = 'pyscript_%s.m' % self.__class__.__name__.split('.')[-1].lower()
        if isdefined(self.inputs.use_mcr) and self.inputs.use_mcr:
            self.mlab.inputs.nodesktop = Undefined
            self.mlab.inputs.nosplash = Undefined
            self.mlab.inputs.single_comp_thread = Undefined
            self.mlab.inputs.uses_mcr = True
            self.mlab.inputs.mfile = True
  
    def _check_mlab_inputs(self):
        if not isdefined(self.inputs.matlab_cmd) and self._matlab_cmd:
            self.inputs.matlab_cmd = self._matlab_cmd
        if not isdefined(self.inputs.paths) and self._paths:
            self.inputs.paths = self._paths
        if not isdefined(self.inputs.use_mcr) and self._use_mcr:
            self.inputs.use_mcr = self._use_mcr

    def _run_interface(self, runtime):
        """Executes the GIFT function using MATLAB."""
        self.mlab.inputs.script = self._make_matlab_command()   	
        results = self.mlab.run()
        runtime.returncode = results.runtime.returncode
        if self.mlab.inputs.uses_mcr:		
            if 'Skipped' in results.runtime.stdout:
                self.raise_exception(runtime)
        runtime.stdout = results.runtime.stdout
        runtime.stderr = results.runtime.stderr
        runtime.merged = results.runtime.merged
        return runtime

    def _list_outputs(self):
        """Determine the expected outputs based on inputs."""
        
        outputs = self._outputs().get()
        return outputs

   
    def _make_matlab_command(self):
        """Generates a mfile to build job structure
       
        Returns
        -------
        mscript : string
            contents of a script called by matlab

        """
        
        raise NotImplementedError

