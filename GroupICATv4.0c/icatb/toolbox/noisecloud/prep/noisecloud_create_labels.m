% Noisecloud: Create image labels

% This is an example of what the .mat file should contain to be used with
% noisecloud.  If you processed your data with FSL, you should use prep/ica

% 1 should be specified as noise, and 0 as not
output.decision = [ 1
                    0
                    1
                    0
                    1 ];
                
output.labels = { 'subject01_IC1'
                  'subject01_IC2'
                  'subject01_IC3'
                  'subject0N_IC1'
                  'subject0N_IC30'};
              
% You will want to save this structure to file, for later use with noisecloud             
save('mylabels.mat')