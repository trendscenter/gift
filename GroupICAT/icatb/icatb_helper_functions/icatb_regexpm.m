function [ cArray ] = icatb_regexpm( varargin )
%+ONELINER:    returns logical t/f for regexp for cell arrays of strings  
% NAME:        regexpm
% FACILITY:    regexp
% SEARCH:      regexp, regexpi, cell array
% LANGUAGE:    matlab
% EXECUTABLE:  yes
% USAGE:       cArray = regexpm( {'top' 'bottom' 'left' 'right'}, {'top' 'left'} );
% URL:         
% AUTHOR:      bettyann chodkowski
% CREATED:     2004-08-01
%
% ARGUMENTS:   
%  str - the string to match; first input argument to regexp
%  expr - the regular expression to match with; second input argument to regexp
%
% RETURNS:     
%  cArray - cell array of logicals
%
% MODIFIED:    2006-03-23
%     2006-03-23 - bettyann chodkowski - chodkowski@kennedykrieger.org
%  added try/catch logic to catch an error that started with the upgrade
%  to R14.  regexp no longer works with cell arrays that are of different
%  lengths.  now an explicit loop is used to make this work.
%
% DESCRIPTION:
%     regexp( cellArray1, cellArray2 ) returns a cell array of start positions
%  and they are awkward to work with.  often i only need to know if there were 
%  matches.  so i run the results of the regexp function thru:
%     ~cellfun( 'isempty', regexp... )
%  to get a logical cell array.
%
% EXAMPLE:
%   <code>
%   >> regexpm( {'top' 'bottom' 'left' 'right'}, {'top' 'left'} )
%   ans =
%        1     0
%        0     0
%        0     1
%        0     0
%   >> find( regexpm( {'top' 'bottom' 'left' 'right'}, {'top' 'left'} ))
%   ans =
%        1
%        7
%   >> any( regexpm( {'top' 'bottom' 'left' 'right'}, {'top' 'left'} )) 
%   ans =
%        1     1
%   </code>
%-

try,
   tmp = regexp( varargin{:} );
catch,
   for ii = [1:length(varargin{1})]
      tmp(ii,:) = regexp( varargin{1}(ii), varargin{2} );
   end;
end;

if ( iscell( tmp ))
   cArray = ~cellfun( 'isempty', tmp );
else
   cArray = ~isempty( tmp );
end;


return;
