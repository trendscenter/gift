function [beta_weights, line_fit] = icatb_constrained_lse(X, data, C)
%% Function to solve constrained least squares problem (y = XB + b) subject to equality
% constraint (CX = 0) using QR decomposition. 
%
% Inputs:
% 1. X - Design matrix
% 2. data - Data
% 3. C - Constraints matrix
% 
% Outputs:
% 1. beta_weights - Beta weights
% 2. line_fit - Line fit
%

[num_rows, num_cols] = size(X);

data = data(:);

dataLength = length(data);

if (num_rows ~= dataLength)
    error('Error:DataLength', 'No. of rows (%d) of design matrix must be the same as the length of the data (%d)\n', num_rows, dataLength);
end

if (size(C, 2) ~= num_cols)
    error('Error:ConstraintSize', 'No. of cols (%d) of constraint matrix must be the same as the no. of cols of design matrix (%d)\n', size(C, 2), num_cols);
end


%% QR decomposition of constraints matrix
[Qc, Rc, Ec] = qr(C');

rank_cmat = rank(Rc);

%% Get null space of the constraints matrix
Qc = Qc(:, rank_cmat + 1:end);

%% Project design matrix to null space
Xproj = X*Qc;

%% Get economy size QR 
[Qx, Rx, Ex] = qr(Xproj, 0);
rank_xmat = rank(Rx);
Qx = Qx(:, 1:rank_xmat);
Rx = Rx(1:rank_xmat, 1:rank_xmat);
Ex = Ex(1:rank_xmat);

%% Fit y to design matrix
line_fit = Qx'*data;       

temp = Rx \ line_fit;

z = zeros(size(Xproj, 2), 1); 

%% Ex contains permutation matrix
z(Ex) = temp;      

%% Get beta coefficients back in full space
beta_weights = Qc * z;             