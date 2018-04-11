% This function prompts the user to select a list of FSL output .ica
% directories and returns the network paths and network row names (the
% identifier of the image).  The subject ID will be extracted from the .ica
% folder name.
%
% Part of the noisecloud package
% vsochat@stanford.edu

function [ network_paths network_row_names ] = noisecloud_read_fsldirs

    selected = 0;
        while selected == 0
            [folders,selected] = spm_select([1 Inf],'dir','Select FSL ICA directories:','','','.ica');
        end

        % Create master list to hold paths and network row names and ICnumbers
        network_paths = [];
        network_row_names = [];

        % Cycle through each folder name, extract all network paths, in order
        for i=1:size(folders,1)

           % Clear the single subject matrices
           [path, fname] = fileparts(folders(i,1:end-1));
           icadir = strcat(path,'\',fname,'.ica');

           % Get list of all zstat images
           zstatholder = dir(strcat(icadir,'\stats\thresh_zstat*.nii.gz'));
           zstatimg = cell(length(zstatholder),1);
           network_row_names_sub = cell(length(zstatholder),1);

           % We are going to make our own file names, 1..length(zstatimg)
           % because matlab's file function makes them out of order!
           for l = 1:length(zstatholder)
               zstatimg{l} = [ icadir '\stats\thresh_zstat' num2str(l) '.nii.gz' ];
               % Now we unzip to create .nii, for use with spm 
               zstatimg(l) = gunzip(zstatimg{l});
               % Note that this will produce the .nii file in the folder
               % This file will be deleted after it is read into matlab.
               network_row_names_sub{l} = [ fname '_IC' num2str(l) ];
           end

           % Now that we have the entire list of images for this subject, add to master list
           network_paths = [ network_paths; zstatimg ];
           network_row_names = [ network_row_names; network_row_names_sub ];
        end
end