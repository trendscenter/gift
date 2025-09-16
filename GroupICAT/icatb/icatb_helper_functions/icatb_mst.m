function [out1 out2 out3] = icatb_mst(A,varargin)
% MST Compute a minimum spanning tree for an undirected graph A.
%
% There are two ways to call MST.
% T = mst(A)
% [i j v] = mst(A)
% The first call returns the minimum spanning tree T of A.
% The second call returns the set of edges in the minimum spanning tree.
% The calls are related by
%    T = sparse(i,j,v,size(A,1), size(A,1));
%    T = T + T';
% The optional algname parameter chooses which algorithm to use to compute
% the minimum spanning tree.  Note that the set of edges returned is not
% symmetric and the final graph must be explicitly symmetrized.
%
% This method works on undirected graphs.
%
% ... = mst(A,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options.
%   options.algname: the minimum spanning tree algorithm
%       ['prim' | {'kruskal'}]
%   options.edge_weight: a double array over the edges with an edge
%       weight for each node, see EDGE_INDEX and EXAMPLES/REWEIGHTED_GRAPHS
%       for information on how to use this option correctly, see Note 1.
%       [{'matrix'} | length(nnz(A)) double vector]
%   options.root: specify the root or starting vertex for the algorithm
%       This option only applies to prim's algorithm.
%       [{'none'} | any vertex number]
%   options.fix_diag: remove any diagonal entries to get correct output
%       from Prim's algorithm [0 | {1}]; beware this option with the
%       edge_weight option too.
%
% Note: the input to this function must be symmetric, so this function
% ignores the 'notrans' default option and never transposes the input.
%
% Note 1: see EXAMPLES/REWEIGHTED_GRAPHS for how to reweight a symmetric
% graph correctly.  There are a few complicated details.
%
% Example:
%    load graphs/clr-24-1.mat
%    mst(A)
%
% See also PRIM_MST, KRUSKAL_MST

% David Gleich
% Copyright, Stanford University, 2006-2008

%% History
%  2006-05-03: Changed to using kruskal as the default following problems
%    with prim due to negative edge weights.
%  2006-05-31: Added full2sparse option
%  2006-06-15: Fixed error with graph symmetric (T+T') instead of max(T,T')
%    found by Mark Cummins
%  2006-11-09: Temporary fix for disconnected graphs and the number of edges
%    in the mst is LESS than n-1.
%  2006-11-10: Added warning for prim with disconnected graphs.
%  2007-04-09: Fixed documentation typos.  (Thanks Chris Maes.)
%  2007-04-09: Fixed bug with 0 weighted graphs.  (Thanks Chris Maes.)
%  2007-04-20: Added edge weight option
%  2007-07-12: Fixed edge_weight documentation
%    Added note about symmetric edge weights
%  2007-12-14: Added rooted option for prim's algorithm
%  2008-10-07: Changed options parsing
%    Addressed issue with incorrect prim output and fixed matrix diagonal
%%

[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end
if trans, end % no trans check

options = struct('algname', 'kruskal', 'edge_weight', 'matrix', ...
    'root', 'none', 'fix_diag', 1);
options = merge_options(options,varargin{:});

fixed_diag= 0;
if options.fix_diag && strcmp(options.algname,'prim') ...
        && any(diag(A)), A = A - diag(diag(A)); fixed_diag= 1; end

% edge_weights is an indicator that is 1 if we are using edge_weights
% passed on the command line or 0 if we are using the matrix.
edge_weights = 0;
edge_weight_opt = 'matrix';

if strcmp(options.edge_weight, 'matrix')
    % do nothing if we are using the matrix weights
else
    edge_weights = 1;
    edge_weight_opt = options.edge_weight;
    if fixed_diag, warning('matlab_bgl:fix_diag',...
            'the diagonal was adjusted, the edge_weight option may be incorrect');
    end
end

if check
    % make sure the matrix is symmetric
    if ~edge_weights
        check_matlab_bgl(A,struct('sym',1,'values',1,...
            'noneg', strcmp(options.algname,'prim')));
    else
        check_matlab_bgl(A,struct());
        
        if strcmp(options.algname,'prim') && any(edge_weights < 0)
            error('matlab_bgl:invalidParameter', ...
                'the edge_weight array must be non-negative');
        end
        
        % check for symmetry
        [i j] = find(A);
        Av = sparse(i,j,edge_weight_opt,size(A,1), size(A,2));
        check_matlab_bgl(Av,struct('sym',1,'nodefault',1));
        
    end
    
    if strcmp(options.algname,'prim')
        if max(components(A)) > 1
            warning('mst:connected', ...
                ['The output from MST using Prim''s algorithm\n' ...
                'on a disconnected graph is only a partial spanning tree.']);
        end
    end
    
end

% old temporary fix for incorrect number of edges
% num_components = max(components(A));

if strcmp(options.root,'none')
    root = 0; % a flag used to denote "no root" to the mex
elseif isa(options.root, 'double')
    root = options.root;
else
    error('matlab_bgl:invalidParameter', ...
        'options.root is not ''none'' or a vertex number.');
end

[i j v] = icatb_mst_mex(A,lower(options.algname),edge_weight_opt,root);

% old temporary fix for disconnected graphs
% if (num_components > 1)
%     i = i(1:end-(num_components-1));
%     j = j(1:end-(num_components-1));
%     v = v(1:end-(num_components-1));
% end

if (nargout == 1 || nargout == 0)
    T = sparse(i,j,v,size(A,1),size(A,1));
    T = T + T';
    out1 = T;
    
    if nnz(T) == 0 && nnz(A) > 0
        warning('mst:empty', ...
            ['MST is empty.  This can occur if you have reweighted\n' ...
            'the matrix with 0 edge weights.  Try the [i j v] = mst(...)\n' ...
            'call instead for that case.']);
    end
else
    out1 = i;
    out2 = j;
    out3 = v;
end;


function check_matlab_bgl(A,options)
% CHECK_MATLAB_BGL Checks the input A for various properties
%
% check_matlab_bgl(A,options) throws an input error if...
%   A is not square
%   if options.values and A is not double valued
%   if options.sym = 1 and A is not symmetric
%   if options.flow_graph = 1 and A is not a flow graph data structure
%   if options.nosparse = 1 do not check if A is sparse
%   if options.nodefault = 1 then do not check default cases
%   if options.nodiag = 1 throw an error if A has any non-zero diagonal values

% David Gleich
% Copyright, Stanford University, 2006-2008

%% History
%  2007-04-20: Added nodefault option
%  2007-07-22: Fixed empty array error for noneg check
%  2008-09-23: Added no diagonal check, misc formatting fixes
%%

if ~isfield(options, 'nodefault') || options.nodefault == 0
    if size(A,1) ~= size(A,2)
        error('matlab_bgl:invalidParameter', 'the matrix A must be square.');
    end
end

if isfield(options, 'values') && options.values == 1
    if ~isa(A,'double')
        error('matlab_bgl:invalidParameter', 'the matrix A must have double values.');
    end
end

if isfield(options, 'noneg') && options.noneg == 1
    v=min(min(A));
    if ~isempty(v) && v < 0
        error('matlab_bgl:invalidParameter', 'the matrix A must have non-negative values.');
    end
end

if isfield(options, 'sym') && options.sym == 1
    if ~isequal(A,A')
        error('matlab_bgl:invalidParameter', 'the matrix A must be symmetric.');
    end
end

if isfield(options, 'nosparse') && options.nosparse == 1
else
    if ~issparse(A)
        error('matlab_bgl:invalidParameter', 'the matrix A must be sparse.  (See set_matlab_bgl_default.)');
    end
end

if isfield(options,'nodiag') && options.nodiag == 1
    if any(diag(A))
        error('matlab_bgl:invalidParameter',...
            'the matrix A must not have any diagonal values')
    end
end



function [trans check full2sparse] = get_matlab_bgl_options(varargin)
%
% Internal private function.
%
% Example:
%    Don't use this function!
%

%% History
%  2008-09-26: Changed to use merge_options instead
%%

doptions = set_matlab_bgl_default();
if nargin>0
    options = merge_options(doptions,varargin{:});
else
    options = doptions;
end

trans = ~options.istrans;
check = ~options.nocheck;
full2sparse = options.full2sparse;

function options = merge_options(default_options,varargin)
% MERGE_OPTIONS Merge a set of default options with options from varargin
% The set of options in varargin can be a list of key,value pairs, or a
% struct with the same information.

% David F. Gleich
% Copyright, Stanford University, 2008

%% History
%  2008-09-25: Initial coding
%%

if ~isempty(varargin) && mod(length(varargin),2) == 0
    options = merge_structs(struct(varargin{:}),default_options);
elseif length(varargin)==1 && isstruct(varargin{1})
    options = merge_structs(varargin{1},default_options);
elseif ~isempty(varargin)
    error('matlag_bgl:optionsParsing',...
        'There were an odd number of key-value pairs of options specified');
else
    options = default_options;
end

function S = merge_structs(A, B)
% MERGE_STRUCTS Merge two structures.
%
% S = merge_structs(A, B) makes the structure S have all the fields from A
% and B.  Conflicts are resolved by using the value in A.
%

%
% merge_structs.m
% David Gleich
%
% Revision 1.00
% 19 Octoboer 2005
%

S = A;

fn = fieldnames(B);

for ii = 1:length(fn)
    if (~isfield(A, fn{ii}))
        S.(fn{ii}) = B.(fn{ii});
    end;
end;

function old_default = set_matlab_bgl_default(varargin)
% SET_MATLAB_BGL_DEFAULT Sets a default option for the Matlab BGL interface
%
% old_default = set_matlab_bgl_default(options) or
% old_default = set_matlab_bgl_default(...) for key-value pair version
% options.istrans: the input matrices are already transposed [{0} | 1]
% options.nocheck: skip the input checking [{0} | 1]
% options.full2sparse: convert full matrices to sparse [{0} | 1]
%
% to get the current set of default options, call
% options = set_matlab_bgl_default()
%
% These options can make the Matlab BGL interface more efficient by
% eliminating the copying operations that occur between Matlab's structures
% and the BGL structures.  However, they are more difficult to use and are
% disabled by default.
%
% Generally, they are best used when you want to perform a large series of
% computations.
%
% Example:
%   % tranpose the matrix initially...
%   At = A'
%   old_options = set_matlab_bgl_default(struct('istrans',1));
%   % perform a bunch of graph work with At...
%   d1 = dfs(At,1); d2 = dfs(At,2); ...
%   % restore the old options
%   set_matlab_bgl_default(old_options);

% David Gleich
% Copyright, Stanford University, 2006-2008

%% History
%%

persistent default_options;
if ~isa(default_options,'struct')
    % initial default options
    default_options = struct('istrans', 0, 'nocheck', 0, 'full2sparse', 0);
end

if nargin == 0
    old_default = default_options;
else
    old_default = default_options;
    default_options = merge_options(default_options,varargin{:});
end

function [ci sizes] = components(A,varargin)
% COMPONENTS Compute the connected components of a graph.
%
% [ci sizes] = components(A) returns the component index vector (ci) and
% the size of each of the connected components (sizes).  The number of
% connected components is max(components(A)).  The algorithm used computes
% the strongly connected components of A, which are the connected
% components of A if A is undirected (i.e. symmetric).
%
% This method works on directed graphs.
% The runtime is O(V+E), the algorithm is just depth first search.
%
% ... = components(A,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options.
%   There are no additional options for this function.
%
% Note: this function does not depend upon the non-zero values of A, but
% only uses the non-zero structure of A.
%
% Example:
%    load('graphs/dfs_example.mat');
%    components(A)
%
% See also DMPERM, BICONNECTED_COMPONENTS

% David Gleich
% Copyright, Stanford University, 2006-2008

%% History
%  2006-04-21: Initial version
%  2006-05-31: Added full2sparse check
%  2006-11-09: Fixed documentation typo.
%  2007-07-08: Code cleanup
%%

[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end
if check, check_matlab_bgl(A,struct()); end
if trans, A = A'; end

[ci sizes] = icatb_components_mex(A);



