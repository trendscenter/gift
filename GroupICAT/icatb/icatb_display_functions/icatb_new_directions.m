function icatb_new_directions(cmd_string);
% load global variable to decide if we want to show directions

% Show the dialog boxes with the appropriate strings
switch lower(cmd_string)
    case 'displaygui-tmap'
        TitleStr = 'Tmap';
        D(1).string = '1. When using Tmaps you need to be careful not to convert the tmap values into z-scores.';
        D(size(D,2)+1).string = '2. Most of the time you do not want to convert tmaps into z-score values.';
        D(size(D,2)+1).string = '3. To leave the t score values unchanged, make sure to convert to z score option is set to no.';
        %str = {D.string};
        %figHandle = icatb_dialogBox('TitleStr', TitleStr, 'textBody', str, 'textType', 'large');


    case 'sortcomponentsgui-about'
        TitleStr = 'Sort Components GUI';

        D(1).string = 'Information On Sorting:';
        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string= 'You have the option of sorting the components using a temporal model or a spatial template.';
        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string = 'Presently there are four types of sorting criteria:';
        D(size(D,2)+1).string = '1. Multiple regression';
        D(size(D,2)+1).string = '2. Correlation';
        D(size(D,2)+1).string = '3. Kurtosis';
        D(size(D,2)+1).string = '4. Maximum Voxel (Only for spatial sorting)';

        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string = 'What do you want to sort?';
        D(size(D,2)+1).string = '1. You can select any subject''s session.';
        D(size(D,2)+1).string = '2. You can select a particuar session of all subjects.';
        D(size(D,2)+1).string = '3. You can select a particuar subject and its sessions.';
        D(size(D,2)+1).string = '4. You can select all the subject''s and sessions.';

        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string = 'Note: ';
        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string = '1. The model timecourses are concatenated to equal the ICA timecourses.';
        D(size(D,2)+1).string = 'You have the option of selecting one model for each subject or one model for all subjects.' ;
        D(size(D,2)+1).string = '2. When more than one regressor is selected for Multiple Linear Regression maintain the same number for all the subjects.';
        D(size(D,2)+1).string = '3. Presently spatial sorting is implemented for any subject any session.';

        %             str = {D.string};
        %             figHandle = icatb_dialogBox('TitleStr', TitleStr, 'textBody', str, 'textType', 'large');

    case 'dispgui-about'
        TitleStr = 'Display GUI';
        D(1).string = 'The display Graphical User Interface(GUI) is used for viewing the results of your group ICA analysis. Display GUI allows you to easily switch viewing methods, viewing sets, and viewing parameters.';
        %str = {D.string};
        %figHandle = icatb_dialogBox('TitleStr', TitleStr, 'textBody', str, 'textType', 'large');


        %Explain how to do single subject analysis
    case 'analysisgui-group-single'
        TitleStr = 'ICA Analysis Directions';
        D(1).string = 'There are three types of ICA analyses; group ICA , single subject single session ICA and single subject multiple session ICA';
        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string =  'A single subject single session analysis is processed the same way that a group analysis is processed except that it only goes through one data reduction step and does not perform other group ICA steps.';
        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string  = 'To run a single subject and single session analysis ->';
        D(size(D,2)+1).string = '1. Select 1 for number of subjects  ';
        D(size(D,2)+1).string = '2. Select 1 for number of sessions  ';
        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string =  'A single subject multiple session analysis is processed the same way that a group analysis is processed except that it does not go through group stats step.';
        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string  = 'To run a single subject and multiple session analysis ->';
        D(size(D,2)+1).string = '1. Select 1 for number of subjects  ';
        D(size(D,2)+1).string = '2. Select more than 1 session  ';
        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string  = 'To run a group analysis ->';
        D(size(D,2)+1).string = '1. Select more than one subject  ';
        D(size(D,2)+1).string = '2. Select any number of sessions ';
        D(size(D,2)+1).string = '';

        %str = {D.string};
        %figHandle = icatb_dialogBox('TitleStr', TitleStr, 'textBody', str, 'textType', 'large');


    case 'analysisgui-datareduction'
        TitleStr = 'Data reduction steps';
        D(1).string = 'Before ICA, the group ICA algorithm uses a data reduction technique called principal component analysis.';
        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string = 'The algorithm works as follows, first a specified number of principal components are extracted from each subject''s functional data. Then the resulting principal components from each subject are concatenated or stacked. This process can be repeated as needed. Once the data has been sufficently reduced the ICA algorithm is run.';
        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string = 'Within the data reduction stage there are a couple of different parameters that you can modify.';
        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string = '1. Number of data reduction steps: The number of times you want to do PCA.';
        D(size(D,2)+1).string = '2. Number of partions to split all the data sets into: The number of partions divided by the number of data sets is the number of data sets in each stack.';
        D(size(D,2)+1).string = '3. Number of principal components to extract';
        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string = 'Note: The number of principal components after the last step of data reduction is the same as the number of independent components you want to extract.';

        %str = {D.string};
        %figHandle = icatb_dialogBox('TitleStr', TitleStr, 'textBody', str, 'textType', 'large');

    case 'analysisgui-calibrate'
        TitleStr = 'Calibrating Timecourses';
        D(1).string = 'After ICA, the timecourses and spatial maps are in arbitary units. Using the original functional data you can scale the components and timecourses to represent percent signal change';

        %str = {D.string};
        %figHandle = icatb_dialogBox('TitleStr', TitleStr, 'textBody', str, 'textType', 'large');


    case 'gift-help'
        TitleStr = 'Group ICA of fMRI Toolbox';
        D(1).string = 'GIFT v1.01b';

        % change later
        %dateStr = date;
        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string = num2str('28-Jul-2004');
        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string = 'GIFT was developed at the Olin Neuropsychiatry Research Center (ONRC) primarily by Vince Calhoun, Eric Egolf and Srinivas Rachakonda with support from NIH R01 EB000840';
        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string = 'Additional contributions were provided by  Baoming Hong and Kent Kiehl at the ONRC, Tulay Adali at the University of Maryland Baltimore County, and Jim Pekar at the FM Kirby Center for Functional Brain Imaging.';
        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string = 'Many additional ICA algorithms were generously contributed by Andrzej Cichocki (also available in the ICALab toolbox at http://www.bsp.brain.riken.jp/ICALAB/ICALABImageProc/).';
        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string = 'GIFT is released under the GNU General Public License';
        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string = 'GIFT uses some functions from the SPM2 library, specifically the spm_vol family of functions and the spm_get function.';
        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string = 'For more information please contact Vince Calhoun at vince.calhoun@yale.edu.';


        %str = {D.string};
        %figHandle = icatb_dialogBox('TitleStr', TitleStr, 'textBody', str, 'textType', 'large');


    case 'convert4dto3d-untested'
        TitleStr = 'Not Well Tested';
        D(1).string = 'The function converts the 4d image file into 3d analyze format. This function is not well tested and still under development.';

        %str = {D.string};
        %figHandle = icatb_dialogBox('TitleStr', TitleStr, 'textBody', str, 'textType', 'large');


    case 'analysis-info'
        TitleStr = 'Analysis Info';
        D(1).string = 'The Analysis Info is used for displaying the information for your group ICA analysis. Analysis Info displays information about parameters, data reduction and output files.';

        %str = {D.string};
        %figHandle = icatb_dialogBox('TitleStr', TitleStr, 'textBody', str, 'textType', 'large');

    case 'runanalysis-info'
        TitleStr = 'Run Analysis';
        D(1).string = 'Purpose: The Run Analysis step consists of:';
        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string = '1. Parameter Initialization';
        D(size(D,2)+1).string = '2. Data Reduction';
        D(size(D,2)+1).string = '3. Calculate ICA';
        D(size(D,2)+1).string = '4. Back Reconstruction';
        D(size(D,2)+1).string = '5. Group Statistics';

        D(size(D,2)+1).string = '';

        D(size(D,2)+1).string = 'Note: Please wait for the analysis to finish as it may take few minutes. ';

        %str = {D.string};
        %figHandle = icatb_dialogBox('TitleStr', TitleStr, 'textBody', str, 'textType', 'large');

    case 'setup-analysis'
        TitleStr = 'Setup ICA Analysis';
        D(1).string = 'Enter the subject information like number of subjects, fMRI data, etc. and then select a suitable ICA algorithm.';

        %str = {D.string};
        %figHandle = icatb_dialogBox('TitleStr', TitleStr, 'textBody', str, 'textType', 'large');

    case 'component-explore'
        TitleStr = 'Component Explorer';
        D(1).string = 'Displays time courses and their image maps in groupings of 1, 4, 9, 16 or 25 per figure. Helps in distinguishing between artifacts and signals of interest.';

        %str = {D.string};
        %figHandle = icatb_dialogBox('TitleStr', TitleStr, 'textBody', str, 'textType', 'large');

    case 'ortho-viewer'
        TitleStr = 'Orthogonal Viewer';
        D(1).string = 'Displays image component overlaid on template Image, BOLD signal for a select voxel, component''s time course and top five component contributors for a given voxel.';

        %str = {D.string};
        %figHandle = icatb_dialogBox('TitleStr', TitleStr, 'textBody', str, 'textType', 'large');

    case 'composite-viewer'
        TitleStr = 'Composite Viewer';
        D(1).string = 'Overlay multiple components of a particular subject''s session or mean for that session.';

        %str = {D.string};
        %figHandle = icatb_dialogBox('TitleStr', TitleStr, 'textBody', str, 'textType', 'large');

    case 'output_files'
        TitleStr = 'Output Files';
        D(1).string = 'All the output files will be prepended with this string.';

    case 'data_files'
        TitleStr = 'Data Files';
        D(1).string = ['No means the functional data needs to be selected. Yes is included as a time saving option if you have already selected files.'];

    case 'mask'
        TitleStr = 'Mask in 3D Analyze Format';
        D(1).string = ['Default mask uses a voxel and compares its intensity value with the mean intensity value of all the voxels in the functional data.', 'Masks must be in the same space as the functional data.'];

    case 'estimate_components'
        TitleStr = 'Dimensionality Estimation';
        D(1).string = ['Components will be automatically estimated from the fMRI data using the MDL criteria.'];

    case 'pc1'
        TitleStr = 'PC1';
        D(1).string = ['Each subject is reduced to the number of principal components selected. For one subject one session PC1 is the number of independent components extracted.'];

    case 'pc2'
        TitleStr = 'PC2';
        D(1).string = ['Each subject is reduced to the number of principal components selected. For data sets less than 4, PC2 is the number of independent components.'];

    case 'pc3'
        TitleStr = 'PC3';
        D(1).string = ['Each subject is reduced to the number of principal components selected. For data sets greater than 10, PC3 is the number of independent components.'];

    case 'z_scores'
        TitleStr = 'Z-Scores';
        D(1).string = ['Scales images and time courses of components to Z-scores after the percent signal change.'];

    case 'adjust_results'
        TitleStr = 'Scaling Results';
        D(1).string = ['After ICA, the timecourses and spatial maps are in arbitary units.'  ...
            'Components can be adjusted by calibrating or converting to z-scores. The explanation of each option is given below:', ...
            'Calibrate: Using the original functional data you can scale the components and timecourses. ', ...
            'Convert to Z-scores: Time courses and spatial maps are converted to Z-scores.'];

    case 'ica_algo'
        TitleStr = 'ICA Algorithms';
        D(1).string = ['ICA options window will popup depending on the algorithm selected. Presently, ICA options are available for Infomax, FastICA and Optimal ICA.'];

    case 'sorting_criteria'
        TitleStr = 'Sorting Criteria';
        D(1).string = ['Components are sorted based on the criteria selected. Multiple Regression, Correlation, Kurtosis are used both for temporal and spatial sorting whereas Maximum Voxel criteria is used only for spatial sorting.'];
    case 'sorting_type'
        TitleStr = 'Sorting Type';
        D(1).string = ['Components can be sorted temporally or spatially. Temporal sorting is a way to compare the model''s time course with the ICA time course whereas spatial sorting classifies the components by comparing the component''s image with the template.'];

    case 'what_to_sort'
        TitleStr = 'What to sort?';
        D(1).string = ['This option is used only for temporal sorting. The options are explained below:'];
        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string = ['1. All Datasets - All data sets are concatenated and correlated with the model time course.'];
        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string = ['2. A Dataset - A particular subject''s session is selected and correlated with the model time course.'];
        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string = ['3. A Session - A particular session for all subjects is selected, time courses are concatenated and correlated with the concatenated model time course.'];
        D(size(D,2)+1).string = '';
        D(size(D,2)+1).string = ['4. A Subject  - A particular subject''s sessions are selected, time courses are concatenated and correlated with the model time course.'];

    otherwise
        TitleStr = 'Error in directions';
        D(1).string = 'Unknown Instructions asked for';

        %str = {D.string};
        %figHandle = icatb_dialogBox('TitleStr', TitleStr, 'textBody', str, 'textType', 'large');
end

str = {D.string};

figHandle = icatb_dialogBox('title', TitleStr, 'textBody', str, 'textType', 'large');

waitfor(figHandle);

%end