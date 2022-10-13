#!/bin/bash
BINDROOT=/data/users2/ceierud/stefan
IN_DIR=$BINDROOT/fmriprep_out
OUT_DIR=$BINDROOT/fmriprep_out_sm
mkdir -p $OUT_DIR
sids=(007 015 018 024 028 044 045 046 051 055 064 067 075 089 093 094 107 112 113 117 120 121 133 134 136 155 166 186 193 198)

for sid in ${sids[@]}; do
	mkdir -p $OUT_DIR/sub-$sid/func/
	fslmaths $IN_DIR/sub-$sid/func/sub-${sid}_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz \
		  -kernel gauss 4.2466452 \
		    -fmean $OUT_DIR/sub-$sid/func/sub-${sid}_task-rest_bold.nii.gz #format filename for BIDS
	mkdir -p $OUT_DIR/sub-$sid/anat/
	cp $IN_DIR/sub-$sid/anat/sub-${sid}_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz $OUT_DIR/sub-$sid/anat/sub-${sid}_T1w.nii.gz
done
