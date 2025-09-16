##Format a raw dataset according to BIDS
# Take it just as an example! Specifics of raw datasets vary greatly across the platforms.
source_path=/data/qneuromark/Data/COBRE/Data_BIDS/Raw_Data/
destination_path=/data/users2/ceierud/stefan/ZN_Neuromark_Raw_BIDS

mkdir -p $destination_path

source_files=($(ls -1 $source_path))
N_subjects=$(($(ls -1 $source_path | wc -l) -1)) #count source_files

#participants.tsv header (required for bids)
echo "subject_id	COBRA_ID" > $destination_path/participants.tsv
#required for bids
echo "{\"Name\": \"ZN Neuromark\", \"BIDS Version\": \"1.0.0rc4\"}" > $destination_path/dataset_description.json
#task desc is required for bids, too
echo "{\"TaskName\": \"DummyTask\", \"RepetitionTime\": 1.0}" > $destination_path/task-rest_bold.json
#README is required for bids
echo "Dummy dataset in BIDS format for a test GIFT-BIDS run" > $destination_path/README 
#Not a requirement for BIDS: Headmotion covariate. We use it later in GIFT dFNC connectivity analysis 
echo "" > headmotion_files.txt
for i in $(seq 1 $N_subjects); do
	subj_id=${source_files[i]} #get subject id from raw dataset
	digit_id=$(printf %03d $i) #3-digit order of subject (001, 002 etc.) 
	echo "${digit_id} ${subj_id}" >> $destination_path/participants.tsv
	mkdir -p $destination_path/sub-$digit_id/func #we call subjects sub-001, sub-002 etc. This is currently needed for GIFT. This may be alleviated in future releases
	cp $source_path/$subj_id/ses_01/func/rest.nii $destination_path/sub-${digit_id}/func/sub-${digit_id}_task-rest_bold.nii # BIDS naming convention of BOLD timeseries
	mkdir -p $destination_path/sub-$digit_id/anat
        cp $source_path/$subj_id/ses_01/anat/T1.nii $destination_path/sub-${digit_id}/anat/sub-${digit_id}_T1w.nii #BIDS naming convention of T1 scans
done


