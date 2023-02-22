function icatb_open_mancovan
%% Open Mancovan Setup GUI
%

param_file = icatb_selectEntry('title', 'Select Mancovan Parameter File', 'typeEntity', 'file', 'typeSelection', 'single', 'filter', '*mancovan.mat');
drawnow;
if (isempty(param_file))
    error('Mancovan parameter file is not selected');
end

icatb_setup_mancovan(param_file);

