function W = zero_diag(W)
%ZERO_DIAG sets the diagonal of a matrix to 0
%   Usage: B = zero_diag(A);
%
%   Input parameters:
%       A   : input matrix
%   Output parameters:
%       B   : output with zero diagonal
%
%   Works also for non-square matrices
%
%   See also: squareform_sp

% code author: Vassilis Kalofolias
% date: March 2016


[m, n] = size(W);

% if not(m==n)
%     warning('non square matrix given!')
% end

n_zeros = min(m, n);

W(1: (m+1): (m+1) * n_zeros) = 0;


