function [parameters] = icatb_get_slice_def(structVol, anatomicalView)
% Purpose: get the transformation matrix slicedef and slices 
% Uses Matt Brett code for display slices
% 
% Input: 
% 1. structVol - Volume of the structural image
% 2. anatomicalView - anatomicalView like axial, sagital or coronal
%
% Output:
% parameters structure with fields 'transformMat', 'slicedef', 'slices',
% 'permuteOrder'.

icatb_defaults;
global USE_DEFAULT_SLICE_RANGE;
global SLICE_RANGE;

useDefault = USE_DEFAULT_SLICE_RANGE;

% Options string
optionsString = {'axial', 'coronal', 'sagittal'};

% get the corresponding index
for ii = 1:length(optionsString)
    index = strmatch(optionsString{ii}, lower(anatomicalView));
    if ~isempty(index)
        orientn = ii;
    end
end

% check the orientn
if isempty(orientn)
    error('Unexpected orientation.');
end

% form transformation matrix
ts = [0 0 0 0 0 0 1 1 1;...
      0 0 0 pi/2 0 0 1 -1 1;...
      0 0 0 pi/2 0 -pi/2 -1 1 1];

% get the transformation matrix
transformMat = icatb_spm_matrix(ts(orientn, :));

% default slice size, slice matrix depends on orientation
% take image sizes from the structural image
D = structVol(1).dim(1:3);

% initialise permute order ([timepoints, sagital, coronal, axial])
permuteOrder = [1 2 3 4];

% permute D1 according to the plane selected
if orientn == 2
    % change permute order
    permuteOrder = [1 2 4 3];
elseif orientn == 3
    % change permute order
    permuteOrder = [1 4 2 3];
end

% transformation matrix
T = transformMat * structVol(1).mat;
vcorners = [1 1 1; D(1) 1 1; 1 D(2) 1; D(1:2) 1; ...
    1 1 D(3); D(1) 1 D(3); 1 D(2:3) ; D(1:3)]';
corners = T * [vcorners; ones(1,8)];
[SC, index] = sort(corners');
vxsz = sqrt(sum(T(1:3,1:3).^2));
slicedef = [SC(1,1) vxsz(1) SC(8,1);SC(1,2) vxsz(2) SC(8,2)];
slices = [SC(1, 3):vxsz(3):SC(8,3)];

% 2006-08-25 - bettyann chodkowski - chodkowski@kennedykrieger.org
% if global variable SLICE_RANGE is defined, then use it
if ( useDefault )
   
   icatb_defaults;
   global SLICE_RANGE;
   if exist( 'SLICE_RANGE', 'var' ) & ~isempty( SLICE_RANGE )
      %sr = strtrim( SLICE_RANGE );
      sr = deblank(SLICE_RANGE);
      sr = regexprep( sr, '\s*', '' );    % remove all whitespace
      sr = regexprep( sr, ':', ' ' );     % replace colons with spaces
      sr = str2num( sr );                 % convert to numeric values
      if ( length(sr) > 2 )
         % SLICE_RANGE includes an increment, ie, a:incr:b
         incr = sr(2);
         sr(2) = [];
      else
         if ( sr(1) < sr(end) ) incr = 1; else incr = -1; end;
      end;
      srList = [ sr(1) : incr : sr(end) ];

      % remove any slices that are not within the valid range of
      % slices defined in 'slices'
      srList(srList < min(slices)) = [];
      srList(srList > max(slices)) = [];

      slices = srList;
   end  
   
end

% output parameter
parameters = struct('transform', transformMat, 'slicedef', slicedef, 'slices', slices, ...
    'permuteOrder', permuteOrder);