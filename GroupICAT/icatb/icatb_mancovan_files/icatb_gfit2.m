function [gf] = icatb_gfit2(t,y,gFitMeasure,options)
% GFIT2 Computes goodness of fit for regression model
%
% USAGE:
%       [gf] = gfit2(t,y)
%       [gf] = gfit2(t,y,gFitMeasure)
%       [gf] = gfit2(t,y,gFitMeasure,options)
%
% INPUT:
%           t:  matrix or vector of target values for regression model
%           y:  matrix or vector of output from regression model.
% gFitMeasure:  a string or cell array of string values representing
%               different form of goodness of fit measure as follows:
%
%               'all' - calculates all the measures below
%               '1' - mean squared error (mse)
%               '2' - normalised mean squared error (nmse)
%               '3' - root mean squared error (rmse)
%               '4' - normalised root mean squared error (nrmse)
%               '5' - mean absolute error (mae)
%               '6' - mean  absolute relative error  (mare)
%               '7' - coefficient of correlation (r)
%               '8' - coefficient of determination (d)
%               '9' - coefficient of efficiency (e)
%               '10' - maximum absolute error
%               '11' - maximum absolute relative error
%
%    options: a string containing other output options, currently the only
%             option is verbose output.
%
%                'v' - verbose output, posts some text output for the
%                      chosen measures to the command line
%
% OUTPUT:
%       gf: vector of goodness of fit values between model output and target for
%       each of the strings in gFitMeasure
%
% EXAMPLES
%
%      gf = gfit2(t,y); for all statistics in list returned as vector
%
%      gf = gfit2(t,y,'3');  for root mean squared error
%
%      gf = gfit2(t,y, {'3'});  for root mean squared error
%
%      gf = gfit2(t,y, {'1' '3' '9'});  for mean squared error, root mean
%            |                          squared error, and coefficient of
%           \|/                         efficiency
%      gf = [mse rmse e]
%
%      gf = gfit2(t,y,'all','v'); for all statistics in list returned as
%                                 vector with information posted to the command
%                                 line on each statistic
%
%      gf = gfit2(t,y, {'1' '3' '9'}, 'v'); for mean squared error, root mean
%                                           squared error, and coefficient of
%                                           efficiency as a vector with
%                                           information on each of these
%                                           also posted to the command line
%
% Modified from gFit 2008/11/07 by Richard Crozier, extended to report multiple
% statistics at once (i.e. can pass in a cell array of strings describing
% desired output), added checking for NaNs, removed superfluous code, added
% max error code and now handles matrices by reshaping them into vectors.
% Also made small improvement to error reporting.
%
% 30 JUN 2009: Added verbose option
% ***********************************************************************
%% INPUT ARGUMENTS CHECK
if (~isempty(which('narginchk')))
    narginchk(2, 4);
else
    error(nargchk(2,4,nargin));
end

% reshape matrices into vectors (order of data is not important)
t = reshape(t,1,[]);
y = reshape(y,1,[]);

if length(t) ~= length(y)
    error('Invalid data size: size of t and y must be same')
end

if nargin > 2
    % check if gFitMeasure is cell string array or just a string
    if ~iscell(gFitMeasure) && ischar(gFitMeasure)
        if strcmp(gFitMeasure,'all')
            % if the string 'all' is passed in, all the stats are
            % required so make the appropriate cell string array
            gFitMeasure = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11'}; % return all measures
        else
            % otherwise convert string to cell string array of size 1
            gFitMeasure = {gFitMeasure};
        end
    else
        % if it is a cell array of strings, check its size
        if size(gFitMeasure,2) == 1
            % if there is only one element check it is not a request
            % for all measures
            if strcmp(char(gFitMeasure),'all')
                % if the string 'all' is passed in, all the stats are
                % required so make the appropriate cell string array
                gFitMeasure = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11'}; % return all measures
            end
        end
    end
else
    % return all measures if only two inputs, nothing will be posted to
    % the command line
    gFitMeasure = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11'};
end

% remove NaNs from the arrays, avoid modifying them if there are no
% NaNs to prevent reallocation of memory
if sum(isnan(t) | isnan(y))
    inds = ~isnan(t) | ~isnan(y);
    t = t(inds);
    y = y(inds);
end

e = t - y; % Calculate the error

gf = ones(1,size(gFitMeasure,2)); % preallocate array

for i = 1:size(gFitMeasure,2)
    
    switch char(gFitMeasure(i))
        
        case '1'                      % mean squared error
            gf(i) = mean(e.^2);        % 0 - perfect match between output and target
            if nargin == 4
                if options == 'v'
                    disp(sprintf(['mean squared error (mse): ' num2str(gf(i))]))
                end
            end
            
        case '2'                      % normalised mean squared error
            gf(i) = mean(e.^2)/var(t); % 0 - perfect match
            if nargin == 4
                if options == 'v'
                    disp(['normalised mean squared error (nmse): ' num2str(gf(i))])
                end
            end
            
        case '3'                      % root mean squared error
            gf(i) = sqrt(mean(e.^2));  % 0 - perfect match
            if nargin == 4
                if options == 'v'
                    disp(['root mean squared error (rmse): ' num2str(gf(i))])
                end
            end
        case '4'                            % normalised root mean squared error
            gf(i) = sqrt(mean(e.^2)/var(t)); % 0 - perfect match
            if nargin == 4
                if options == 'v'
                    disp(['normalised root mean squared error (nrmse): ' num2str(gf(i))])
                end
            end
        case '5'                      % mean absolute error
            gf(i) = mean(abs(e));       % 0 - perfect match
            if nargin == 4
                if options == 'v'
                    disp(['mean absolute error (mae): ' num2str(gf(5))])
                end
            end
        case '6'                      % mean absolute relative error
            gf(i) = mean((abs(e./t)));  % 0 - perfect match
            if nargin == 4
                if options == 'v'
                    disp(['mean  absolute relative error (mare): ' num2str(gf(6))])
                end
            end
        case '7'                      % coefficient of correlation
            cf = corrcoef(t,y);      % 1 - perfect match
            gf(i) = cf(1,2);
            if nargin == 4
                if options == 'v'
                    disp(['coefficient of correlation (r): ' num2str(gf(7))])
                end
            end
        case '8'                      % coefficient of determination (r-squared)
            cf = corrcoef(t,y);
            gf(i) = cf(1,2);
            gf(i) = gf(i)^2;               % 1 - perfect match
            if nargin == 4
                if options == 'v'
                    disp(['coefficient of determination (r-squared): ' num2str(gf(8))])
                end
            end
        case '9'                                       % coefficient of efficiency
            gf(i) = 1 - sum(e.^2)/sum((t - mean(t)).^2); % 1 - perfect match
            if nargin == 4
                if options == 'v'
                    disp(['coefficient of efficiency (e): ' num2str(gf(9))])
                end
            end
        case '10'
            gf(i) = max(abs(e));       % maximum absolute error
            if nargin == 4
                if options == 'v'
                    disp(['maximum absolute error: ' num2str(gf(10))])
                end
            end
        case '11'
            gf(i) = max(abs(e./t));    %  maximum absolute relative error
            if nargin == 4
                if options == 'v'
                    disp(['maximum absolute relative error: ' num2str(gf(11))])
                end
            end
        otherwise
            error('Invalid goodness of fit measure in gFitMeasure(%d):\nIt must be one of the strings {1 2 3 4 5 6 7 8 9 10 11}, but actually contained ''%s''',i,char(gFitMeasure(i)))
    end
end
end
%**************************************************************************










