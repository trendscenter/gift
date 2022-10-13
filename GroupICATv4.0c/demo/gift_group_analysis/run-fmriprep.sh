module load singularity/3.10.2 #cluster-specific. make sure singularity is available on command "singularity"
BINDROOT=/data/users2/ceierud/stefan/ # root directory
BINDP=/data #binding point for root directory
PARTICIPANTS_LABELS=(075)
IN_DIR=$BINDROOT/ZN_Neuromark_Raw_BIDS
OUT_DIR=$BINDROOT/fmriprep_out
WORKDIR=$BINDROOT/fmriprep_work
#freesurfer license
export SINGULARITYENV_FS_LICENSE=/data/users2/ceierud/stefan/container_run_data/license.txt
FMRIPREP_CONTAINER=singularity_containers/fmriprep_latest.sif

#run fMRIprep. it defaults to using all cores and all available memory, to my best knowledge
singularity run -B $BINDP:$BINDP --cleanenv $FMRIPREP_CONTAINER $IN_DIR $OUT_DIR \
	participant --participant-label ${PARTICIPANTS_LABELS[@]} \
	--omp-nthreads=16 --work-dir $WORKDIR
