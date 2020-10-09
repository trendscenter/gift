function g = icatb_imgrescale(data, newsize)
% IMGRESCALE   Rescales raw image data using bilinear interpolation.
%    G = IMGRESCALE(DATA, NEWSIZE) Rescales image data, provided as
%    a MxN matrix or as a MxNx3 matrix (for RGB image data) to a
%    new size NEWSIZE = [W, H] of width W and height H. Both double
%    arrays and uint8 arrays are supported.

% Author : Andreas Klimke, Universität Stuttgart
% Version: 1.0
% Date   : June 17, 2003
	
% get the size of the original image
oldsize = size(data);

% check the data type of the array
oldtype = [];
if ~isa(data, 'double')
	if isa(data, 'uint8')
		oldtype = class(data);
		data = double(data);
	else
		error('Only double or uint8 data types supported.');
	end
end

% check if three color planes are provided or not
if length(oldsize)>2 
	if oldsize(3) == 3
		% RGB data present, otherwise treat as gray scale data.
		RGB = 1;
	else
		error('Format not supported.');
	end
else
	RGB = 0;
end

% compute scaling factors
factor = (oldsize(1:2)-1)./(newsize-1);

% create new grid
u = 0:newsize(1)-1;
v = 0:newsize(2)-1;
[U, V] = ndgrid(u, v);

% transform new grid vectors to old grid size
u = u.*factor(1) + 1;
v = v.*factor(2) + 1;

% Compute the location of each new point relative to one nearest
% neighbor of the original image
U = U.*factor(1); U = U - fix(U);
V = V.*factor(2); V = V - fix(V);

% Perform interpolation element by element
if ~RGB
	g = (V-1).*((U-1).*data(floor(u), floor(v)) - ...
							U.*data(ceil(u), floor(v))) - ...
			V.*((U-1).*data(floor(u), ceil(v)) - ...
					U.*data(ceil(u), ceil(v)));
else
	% replicate grid to all three dimensions
	U = repmat(U, [1 1 3]);
	V = repmat(V, [1 1 3]);
	g = (V-1).*((U-1).*data(floor(u), floor(v), :) - ...
							U.*data(ceil(u), floor(v), :)) - ...
			V.*((U-1).*data(floor(u), ceil(v), :) - ...
					U.*data(ceil(u), ceil(v), :));
end	

% reconvert to old data format, if necessary
if ~isempty(oldtype)
	switch oldtype
	 case 'uint8'
		g = uint8(g);
	end
end
