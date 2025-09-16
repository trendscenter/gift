module load singularity/3.10.2 #cluster-specific. make sure singularity is available on command "singularity"
BINDROOT=/data/users2/ceierud/stefan/ # root directory
BINDP=/data #binding point for root directory
WORKDIR=$BINDROOT/fmriprep_work_sm
SMOOTH_SCRIPT=$BINDROOT/smooth_subjects.sh
mkdir -p $WORKDIR

#freesurfer license
export SINGULARITYENV_FS_LICENSE=/data/users2/ceierud/stefan/container_run_data/license.txt
FMRIPREP_CONTAINER=singularity_containers/fmriprep_latest.sif

singularity exec -B $BINDP:$BINDP --cleanenv $FMRIPREP_CONTAINER $SMOOTH_SCRIPT

#add anatomical & dataset description
cd ZN_Neuromark_Raw_BIDS
cp task-rest_bold.json dataset_description.json participants.tsv README ../fmriprep_out_sm/
cd ..
