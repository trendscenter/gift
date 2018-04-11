function ICA_Options = icatb_sbica_options(ICA_Options, dewhiteM)
% Add some more options to use semi-blind ICA

% Time dimension
tdim = size(dewhiteM, 1);

whiteM = pinv(dewhiteM);

% Number of Independent Comp.
numOfIC = size(dewhiteM, 2);

%tar and nov are SPM time courses
RA = randn(tdim, numOfIC)/.05;

ind = strmatch('TC', ICA_Options(1:2:end), 'exact');

if isempty(ind)
    error('Time course constraints are not selected.');
end

TC = ICA_Options{2*ind(1)};

% Make the timecourses to match the tdim (number of time points) by ICdim
% (number of components)
if size(TC, 1) ~= tdim
    TC = TC';
end

% Average of the time course constraints
RA(:, 1) = sum(TC, 2) / size(TC, 2);


all_opt_names = ICA_Options(1:2:end);
[ddd, inds] = unique(all_opt_names);
good_inds = (1:length(ICA_Options));
good_inds = sort(good_inds([2.*inds - 1, 2.*inds]));
ICA_Options = ICA_Options(good_inds);


% get the number of ICA options
%numOptions = length(ICA_Options);

ICA_Options = addOpts(ICA_Options, 'weights', inv(whiteM*RA)/10000000);

chkInds = strmatch('prefs', lower(ICA_Options(1:2:end)), 'exact');

if (isempty(chkInds))
    % Pass dewhiteM, prefs
    prefs = [0.5, zeros(1, numOfIC - 1)];
    ICA_Options{end + 1} = 'prefs';
    ICA_Options{end + 1} = prefs;
else
    prefs = ICA_Options{2*chkInds};
    prefs = [prefs(1), zeros(1, numOfIC - 1)];
    ICA_Options{2*chkInds} = prefs;
end

ICA_Options = addOpts(ICA_Options, 'dewhiteM', dewhiteM);
ICA_Options = addOpts(ICA_Options, 'whiteM', whiteM);


function ICA_Options = addOpts(ICA_Options, str, val)

chkInds = strmatch(str, lower(ICA_Options(1:2:end)), 'exact');

if (isempty(chkInds))
    ICA_Options{end + 1} = str;
    ICA_Options{end + 1} = val;
else
    ICA_Options{2*chkInds} = val;
end