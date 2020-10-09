function [sR] = icatb_icassoEst(mode, X, M, varargin)
%function sR=icassoEst(mode,X,M,['FastICAparamName1',value1,'FastICAparamName2',value2,...])
%
%PURPOSE
%
%To compute randomized ICA estimates M times from data X. Output of
%this function (sR) is called 'Icasso result structure' (see
%icassoStruct). sR keeps the on all the methods, parameters, and
%results in the Icasso procedure.
%
%EXAMPLES OF BASIC USAGE
%
%   sR=icassoEst('randinit', X, 30);
%
%estimates ICA for 30 times on data matrix X using Icasso
%default parameters for FastICA: symmetrical approach, kurtosis as
%contrast function. In maximum 100 iterations are used for
%estimating ICA in each round. Randomizes only initial conditions.
%
%   sR=icassoEst('both', X, 30, 'g', 'tanh', 'approach', 'defl');
%
%estimates ICA for 15 times on data matrix X using 'tanh' as the
%contrast function and the deflatory approach in FastICA. Applies
%both bootstrapping the data and randomizing initial conditions.
%
%INPUT
%
% mode (string) 'randinit' | 'bootstrap | 'both'
% X    (dxN matrix) data where d=dimension, N=number of vectors
% M    (scalar) number of randomizations (estimation cycles)
%
%Optional input arguments are given as argument identifier - value
%pairs: 'identifier1', value1, 'identifier2', value2,...
%(case insensitive)
%
%FastICA parameters apply here (see function fastica)
%Default: 'approach', 'symm', 'g', 'pow3', 'maxNumIterations', 100
%
%OUTPUTS
%
% sR (struct) Icasso result data structure
%
%DETAILS
%
%Meaning of different choices for input arg. 'mode'
% 'randinit': different random initial condition each time.
% 'bootstrap': the same initial cond. each time, but data is
%   bootstrapped. The initial condition can be explicitly
%   specified using FastICA parameter 'initGuess'.
% 'both': use both data bootstrapping and randomization of
%    initial condition.
%
%FASTICA PARAMETERS See function 'fastica' in FastICA toolbox for
%more information. Note that the following FastICA parameters
%cannot be used:
%
% In all modes ('randinit','bootstrap','both'):
%   using 'interactivePCA','sampleSize', 'displayMode', 'displayInterval',
%   and 'only' are not allowed for obvious reasons. In addition,
% in modes 'randinit' and 'both':
%   using 'initGuess' is not allowed since initial guess is
%   randomized, and
% in modes 'bootstrap' and 'both':
%   using 'whiteMat', 'dewhiteMat', and 'whiteSig' are not allowed
%   since they need to be computed for each bootstrap sample
%   individually.
%
%ESTIMATE INDEXING CONVENTION: when function icassoEst is run
%each estimate gets a unique, integer label in order of
%appearance. The same order and indexing is used throughout the
%Icasso software. In many functions, one can pick a subset of
%estimates sR by giving vector whose elements refers to this unique
%label.
%
%SEE ALSO
% icasso
% fastica
% icassoStruct
% icassoExp
% icassoGet
% icassoShow
% icassoResult
%
%When icassoEst is accomplished, use icassoExp to obtain clustering
%results and to store them in sR After this, the results can be
%examined visually using icassoShow. Results and other information
%an be finally retrieved also by functions icassoResult and icassoGet.

%COPYRIGHT NOTICE
%This function is a part of Icasso software library
%Copyright (C) 2003-2005 Johan Himberg
%
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% ver 1.1 johan 210704

% Set the Icasso struct


% varargin contains the following whiteM, dewhiteM and whitesig
% algorithm name, icaOptions

if (isa(X, 'single'))
    X = double(X);
end

sR = icassoStruct(X);

%mode = 'randinit';

sR.mode = mode;

for ii = 1:2:length(varargin)
    if strcmpi(varargin{ii}, 'whitem')
        White = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'dewhitem')
        deWhite = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'whitesig')
        w = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'algoindex')
        algoIndex = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'icaoptions')
        icaOptions = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'numofpc')
        numOfPC = varargin{ii + 1};
    end
end

if (isa(White, 'single'))
    White = double(White);
end

if (isa(deWhite, 'single'))
    deWhite = double(deWhite);
end

% store whitening for original data:
sR.whiteningMatrix = White;
sR.dewhiteningMatrix = deWhite;

icaAlgo = icatb_icaAlgorithm;
if ischar(algoIndex)
    algoIndex = strmatch(lower(algoIndex), lower(icaAlgo), 'exact');
end

algorithmName = deblank(icaAlgo(algoIndex, :));


%% Compute N times FastICA
k=0; index=[];
for i=1:M,
    %clc;
    fprintf('\n\n%s\n\n',['Randomization using ', algorithmName, ': Round ' num2str(i) '/' ...
        num2str(M)]);

    switch mode
        case 'randinit'
            % data is fixed;
            X_=X;
        case {'bootstrap','both'}
            % Bootstrap and compute whitening for _bootstrapped_ data
            X_=bootstrap(X);
            %[w,White,deWhite]=fastica(X_,'only','white',fasticaoptions{:});
            %%%%%% Calculate PCA and Whitening matrix %%%%%
            % PCA
            [V, Lambda] = icatb_v_pca(X_, 1, numOfPC, 0, 'transpose', 'no');

            % whiten matrix
            [w, White, deWhite] = icatb_v_whiten(X_, V, Lambda, 'transpose');

            %%%%% End for calculating ICA and whitening matrix %%%%
        otherwise
            error('Internal error?!');
    end

    % Estimate FastICA set displayMode off
    %   [dummy,A_,W_]=fastica(X_,fasticaoptions{:},...
    % 			'whiteMat',White,'dewhiteMat',deWhite,'whiteSig',w,...
    % 			'sampleSize',1,'displayMode','off');

    [icaAlgo, W_, A_] = icatb_icaAlgorithm(algoIndex, w, icaOptions);

    A_ = deWhite*pinv(W_);
    W_ = pinv(A_);

    % Store results if any
    n=size(A_,2);
    if n>0,
        k=k+1;
        sR.index(end+1:end+n,:)=[repmat(k,n,1), [1:n]'];
        sR.A{k}=A_; sR.W{k}=W_;
    end
end


function X=bootstrap(X)

N=size(X,2);
index=round(rand(N,1)*N+.5);
X=X(:,index);