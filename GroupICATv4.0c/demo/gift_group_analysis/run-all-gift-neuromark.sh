# Run GIFT-BIDS on all subjects of a study using a Singularity container (converted from a Docker conainer)
module load singularity/3.10.2
ROOT=/data/users2/ceierud/stefan
TR_PWD=$(pwd) # working directory of the script and container
DATASET_NAME=fmriprep_out_sm
TIMESTAMP=$(date +"%y%m%d%H%M")
LOG=gift_subj001-203_$TIMESTAMP.log
OUT_DIR=MF_sm_GIFT
sids=(007 015 018 024 028 044 045 046 051 055 064 067 075 089 093 094 107 112 113 117 120 121 133 134 136 155 166 186 193 198)
mkdir -p $OUT_DIR

mkdir -p $PWD/tmp1
mkdir -p $PWD/tmp2
echo "Starting run: $(date +"%H:%M:%S")"
# Run GIFT-BIDS container in Singularity
singularity run --bind $TR_PWD/tmp1:/tmp --bind $TR_PWD/tmp2:/var/tmp --bind ${TR_PWD}/${DATASET_NAME}:/data --bind ${TR_PWD}/${OUT_DIR}:/output --bind $TR_PWD/container_run_data/cfg:/cfg $ROOT/singularity_containers/trends_gift-bids-2022-05-23-9c12f738c4dd.img /data /output participant --participant_label ${sids[@]} --config /cfg/config_spatial_ica_bids.m 1>$LOG 2>&1
echo "Finished: $(date +"%H:%M:%S")"

