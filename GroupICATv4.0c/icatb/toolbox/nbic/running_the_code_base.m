%% An example runner for the code base. The follwoin steps will run the biclsutering analysis
%% Step 1:  Load the raw data matrix and specify the user arguments   
clear
clc
load ('Input_data.mat')
% Specify the user preferences for the biclusters. Standard values set 
x_size = 50; % minimum number of samples 
y_size = 3;  % minimum number of features/variables 
tolerance = 20; % Percentage of allowed overlap between two biclusters 
reps = 5; % repititions
% Create Feature IDs. For more flexibility in skipping some features from the analysis 
variable_Ids = 1:size(raw_data_matrix,2);
%% Step 2: Run the main subroutine 'NBiC' with relevant arguments specified above. It will call all the require subroutines 
% Returns a set of biclusters.
% Each bicluster has three fields. 
% subs: number of subjects inlcuded
% comps: number of components/features inlcuded
% freq: stability across the repititions. 
% 'biclusters' is the desired resul - a set of biclusters. You need to store it in the memory  
biclusters = NBiC(x_size, y_size, tolerance, variable_Ids, reps, raw_data_matrix);

% save ('biclusters_new.mat','biclusters')

%%