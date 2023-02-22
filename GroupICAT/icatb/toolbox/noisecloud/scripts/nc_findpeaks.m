function [pks,locs] = nc_findpeaks(X,varargin)

[X,Ph,Pd,Th,Np,Str,infIdx] = parse_inputs(X,varargin{:});
[pks,locs] = getPeaksAboveMinPeakHeight(X,Ph);
[pks,locs] = removePeaksBelowThreshold(X,pks,locs,Th,infIdx);
[pks,locs] = removePeaksSeparatedByLessThanMinPeakDistance(pks,locs,Pd);
[pks,locs] = orderPeaks(pks,locs,Str);
[pks,locs] = keepAtMostNpPeaks(pks,locs,Np);

%--------------------------------------------------------------------------
function [X,Ph,Pd,Th,Np,Str,infIdx] = parse_inputs(X,varargin)

% % Validate input signal
% validateattributes(X,{'numeric'},{'nonempty','real','vector'},...
%     'findpeaks','X');
% try
%     % Check the input data type. Single precision is not supported.
%     chkinputdatatype(X);
% catch ME
%     throwAsCaller(ME);
% end
M = numel(X);
if (M < 3)
    error('empty dataset');
end


Ph  = -Inf;
Pd  = [];
Th  = 0;
Np  = [];
Str = 'none';



% Set default values for MinPeakDistance and NPeaks
if(isempty(Pd)), Pd = 1; end
if(isempty(Np)), Np = M; end

if(Pd >= M)
    error(message('signal:findpeaks:largeMinPeakDistance', 'MinPeakDistance', 'MinPeakDistance', num2str( M )));
end

% Replace Inf by realmax because the diff of two Infs is not a number
infIdx = isinf(X);
if any(infIdx),
    X(infIdx) = sign(X(infIdx))*realmax;
end
infIdx = infIdx & X>0; % Keep only track of +Inf

%--------------------------------------------------------------------------
function [pks,locs] = getPeaksAboveMinPeakHeight(X,Ph)

pks = [];
locs = [];

if all(isnan(X)),
    return,
end

Indx = find(X > Ph);
if(isempty(Indx))
    warning(message('signal:findpeaks:largeMinPeakHeight', 'MinPeakHeight', 'MinPeakHeight'));
    return
end
    
% Peaks cannot be easily solved by comparing the sample values. Instead, we
% use first order difference information to identify the peak. A peak
% happens when the trend change from upward to downward, i.e., a peak is
% where the difference changed from a streak of positives and zeros to
% negative. This means that for flat peak we'll keep only the rising
% edge.
trend = sign(diff(X));
idx = find(trend==0); % Find flats
N = length(trend);
for i=length(idx):-1:1,
    % Back-propagate trend for flats
    if trend(min(idx(i)+1,N))>=0,
        trend(idx(i)) = 1; 
    else
        trend(idx(i)) = -1; % Flat peak
    end
end
        
idx  = find(diff(trend)==-2)+1;  % Get all the peaks
locs = intersect(Indx,idx);      % Keep peaks above MinPeakHeight
pks  = X(locs);

%--------------------------------------------------------------------------
function [pks,locs] = removePeaksBelowThreshold(X,pks,locs,Th,infIdx)

idelete = [];
for i = 1:length(pks),
    delta = min(pks(i)-X(locs(i)-1),pks(i)-X(locs(i)+1));
    if delta<Th,
        idelete = [idelete i]; %#ok<AGROW>
    end
end
if ~isempty(idelete),
    locs(idelete) = [];
end

X(infIdx) = Inf;                 % Restore +Inf
locs = union(locs,find(infIdx)); % Make sure we find peaks like [realmax Inf realmax]
pks  = X(locs);

%--------------------------------------------------------------------------
function [pks,locs] = removePeaksSeparatedByLessThanMinPeakDistance(pks,locs,Pd)
% Start with the larger peaks to make sure we don't accidentally keep a
% small peak and remove a large peak in its neighborhood. 

if isempty(pks) || Pd==1,
    return
end

% Order peaks from large to small
[pks, idx] = sort(pks,'descend');
locs = locs(idx);

idelete = ones(size(locs))<0;
for i = 1:length(locs),
    if ~idelete(i),
        % If the peak is not in the neighborhood of a larger peak, find
        % secondary peaks to eliminate.
        idelete = idelete | (locs>=locs(i)-Pd)&(locs<=locs(i)+Pd); 
        idelete(i) = 0; % Keep current peak
    end
end
pks(idelete) = [];
locs(idelete) = [];

%--------------------------------------------------------------------------
function [pks,locs] = orderPeaks(pks,locs,Str)

if isempty(pks), return; end

if strcmp(Str,'none')
    [locs,idx] = sort(locs);
    pks = pks(idx);
else
    [pks,s]  = sort(pks,Str);
    locs = locs(s);
end

%--------------------------------------------------------------------------
function [pks,locs] = keepAtMostNpPeaks(pks,locs,Np)

if length(pks)>Np,
    locs = locs(1:Np);
    pks  = pks(1:Np);
end

% [EOF]
