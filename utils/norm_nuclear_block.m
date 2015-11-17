function [n] = norm_nuclear_block(X, ind_rows, ind_cols, param)
%NORM_NUCLEAR_BLOCK Blocwise nuclear norm
%   Usage:  n = norm_nuclear_block(X, ind_rows, ind_cols)
%
%   Input parameters:
%           X     : Input matrix
%           ind_r : Vector partitioning the rows of X in groups
%           ind_c : Vector partitioning the columns in groups 
%           param : Optional paramter
%   Output parameters:
%           n     : Norm
%           info  : Structure summarizing information at convergence
%
%   `prox_NuclearNorm_Block(x, gamma, param)` solves:
%
%   .. n = sum_{i,j} W(i,j)||Z(i,j)||_*
%
%   .. math::  n =  \sum_{i,j} ||Z_{i,j}||_*
%   
%   where Z(i,j) is the i,j-th block indicated by the indices ind_r == i,
%   ind_c == j and w(i,j) is an optional weight for the block
%
%   EXAMPLE: ind_r [1 2 2 3 3 1] means that the first block
%   contains the first and last rows of x
%
%   * *param.single* : single precision (1) or not (0)? (default: single
%     only if input is single precision);
%
%   See also:  prox_nuclearnorm norm_nuclear 
%

% code author: Nathanael Perraudin,Vassilis Kalofolias
% date: November 2015



if nargin < 4, param = struct; end
if ~isfield(param, 'single'), param.single = isa(X, 'single'); end

Nr = max(ind_rows);
Nc = max(ind_cols);





n = 0;

for i_r = 1:Nr
    for i_c = 1:Nc
        % compute the SVD of each small block      
        if param.single
            sigma = svd(single(full(X(ind_rows == i_r, ind_cols == i_c))), 'econ');       % good for small, dense matrices!!
        else
            sigma = svd(full(X(ind_rows == i_r, ind_cols == i_c)), 'econ');       % good for small, dense matrices!!
        end
        % column vector  
%         sigma = diag(S);                     
        nuclearNorm =  sum(abs(sigma));

        n = n + nuclearNorm;
        
        % keep statistics separately for each block if asked
    end
end




end


