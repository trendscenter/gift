Script you may use to run Neuromark 2.2 together with GIFT v 4.0.6.31 or later.

Steps to run
1. Download GIFT version 4.0.6.31 or later from https://github.com/trendscenter/gift
2. Open matlab
3. In matlab command line type: addpath(genpath('c:\where\you\downloaded\gift\GroupICAT\icatb'))
4. In matlab command line type: copyfile(which('input_neuromark_fmri_2.2.m'),'c:\where\my\project\is')
5. In matlab command line type: edit('c:\where\my\project\is\input_neuromark_fmri_2.2.m')
5. In the input file (now in edit mode), change line for outputDir so: "outputDir = 'c:\where\my\project\is\output';"
6. In same file, add your subjects in input_data_file_patterns, with example for two subjects: "input_data_file_patterns = {'c:\where\my\project\is\input\subject1.nii',    'c:\where\my\project\is\input\subject2.nii'};"
7. Of course you can change this to the amount of subjects you want to run
8. Make sure to save the file you have edited in steps above
9. In matlab command line type: icatb_batch_file_run('c:\where\my\project\is\input_neuromark_fmri_2.2.m');
11. You should get results and reports
