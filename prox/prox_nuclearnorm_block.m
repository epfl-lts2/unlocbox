function [sol, info] = prox_nuclearnorm_block(X, gamma, ind_rows, ind_cols, param)
%PROX_NUCLEARNORM_BLOCK Proximal operator of nuclear norms of blocks
%   Usage:  sol = prox_nuclearnorm_block(x, gamma, ind_r, ind_c)
%           sol = prox_nuclearnorm_block(x, gamma, ind_r, ind_c, param)
%           [sol, info] = prox_nuclearnorm_block(...)
%
%   Input parameters:
%           X     : Input matrix
%           gamma : Regularization parameter
%           ind_r : Vector partitioning the rows of X in groups
%                   EXAMPLE: ind_r [1 2 2 3 3 1] means that the first block
%                   contains the first and last rows of x
%           ind_c : Vector partitioning the columns in groups (same as
%                   ind_r)
%           param : Structure of optional parameters
%   Output parameters:
%           sol   : Solution
%           info  : Structure summarizing information at convergence
%
%   `prox_NuclearNorm_Block(x, gamma, param)` solves:
%
%   .. sol = argmin_{Z} 0.5*||X - Z||_F^2 + sum_{i,j} gamma * W(i,j)||Z(i,j)||_*
%
%   .. math::  sol = arg\min_{Z} \frac{1}{2} \|X - Z\|_F^2 + \sum_{i,j} \gamma W_{i,j}  ||Z_{i,j}||_*
%   
%   where Z(i,j) is the i,j-th block indicated by the indices ind_r == i,
%   ind_c == j and w(i,j) is an optional weight for the block
%
%   param is a Matlab structure containing the following fields:
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence, 2 print info
%     for each block (default: 1)
%
%   * *param.single* : single precision (1) or not (0)? (default: single
%     only if input is single precision);
%   
%   * *param.compute_stat* : if true, the statistics nz_blocks, rank_block,
%     norm_block will be returned as fields of the struct info.
%
%   * *param.W* : weight for the term of each block in form of an array.
%
%   info is a Matlab structure containing the following fields:
%
%   * *info.algo* : Algorithm used
%
%   * *info.iter* : Number of iteration
%
%   * *info.time* : Time of execution of the function in sec.
%
%   * *info.final_eval* : Final evaluation of the sum of nuclear norms
%
%   * *info.crit* : Stopping criterion used 
%
%   * *info.rank* : Rank of the final solution (-1 means the rank was not
%     computed) 
%   
%   * *info.nz_blocks* :    total number of zero blocks
%
%   * *info.rank_block* :   array containing the rank of each block
%
%   * *info.norm_block* :   array containing the nuclear norm of each block
%
%
%   See also:  prox_nuclearnorm prox_l1 proj_b1 prox_tv

%
% code author: Vassilis Kalofolias
% date: Feb 2015


if nargin < 5, param = struct; end

if ~isfield(param, 'single'), param.single = isa(X, 'single'); end
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'compute_stat'), 
    compute_stat = false;
else
    compute_stat = param.compute_stat;
end
if ~isfield(param, 'W')
    weights_given = false;
else
    weights_given = true;
end


t1 = tic;

Nr = max(ind_rows);
Nc = max(ind_cols);

% Test of gamma
if gamma == 0
    sol = X;
    % set the information structure returned
    info.algo = mfilename;
    info.iter = 0;
    info.final_eval = 0;
    info.crit = '--';
    info.time = toc(t1);
    info.rank_sum = -1;     %rank(full(sol));
    if compute_stat
        info.nz_blocks = -1;    % number of zero blocks
        info.rank_block = -1 * ones(Nr, Nc);  % rank of each block
        info.norm_block = -1 * ones(Nr, Nc);  % nuclear norm of each block
    end
    return
end




% Useful functions
soft = @(z, T) sign(z) .* max(abs(z) - T, 0);

[m, n] = size(X);

if length(ind_rows) ~= m || length(ind_cols) ~= n
    error('Lengths of indices of rows and columns have to span the whole set of rows and columns of input matrix');
end

sol = X;
if compute_stat
    rank_block = zeros(Nr, Nc);
    norm_block = zeros(Nr, Nc);  % rank of each block
end

nuclearNorm_sum = 0;
rank_sum = 0;
nz_blocks = 0;
weight_ij = 1;

for i_r = 1:Nr
    for i_c = 1:Nc
        % compute the SVD of each small block
        if weights_given
            weight_ij = param.W(i_r, i_c);
        end
        
        if param.single
            [U, S, V] = svd(single(full(X(ind_rows == i_r, ind_cols == i_c))), 'econ');       % good for small, dense matrices!!
        else
            [U, S, V] = svd(full(X(ind_rows == i_r, ind_cols == i_c)), 'econ');       % good for small, dense matrices!!
        end
       
        % Shrink:
        sigma = diag(S);                % column vector
        sigma = soft(sigma, gamma * weight_ij);     % modified singular values!
        r = sum(sigma > 0);             % rank of solution
        
        % Reconstruct X with new singular values
        nuclearNorm = weight_ij * sum(sigma(1:r));
        % sol = Ur Sr Vr', where Ur = U(:, 1:r), Sr = diag(sigma), Vr = Vr(:, 1:r)
        sol(ind_rows == i_r, ind_cols == i_c) = U(:, 1:r) * bsxfun(@times, sigma(1:r), V(:, 1:r).');
        
        if param.verbose >= 2
            fprintf('prox_nuclearnorm_block: block (%3i, %3i):  rank = %3i,  |X|_* = %5e \n', i_r, i_c, r, nuclearNorm);
        end
        
        nuclearNorm_sum = nuclearNorm_sum + nuclearNorm;
        rank_sum = rank_sum + r;
        nz_blocks = nz_blocks + (r == 0);
        
        % keep statistics separately for each block if asked
        if compute_stat
            rank_block(i_r, i_c) = r;
            norm_block(i_r, i_c) = nuclearNorm;
        end
    end
end


if param.verbose > 0
    fprintf('prox_nuclearnorm_block: sum|X|_* = %5e, sum_rank = %3i,  %3d zero blocks:  \n', nuclearNorm_sum, rank_sum, nz_blocks);
end

% set the information structure returned
info.algo = mfilename;
info.iter = 0;
info.final_eval = nuclearNorm_sum;
info.crit = '--';
info.rank_sum = rank_sum;     %rank(full(sol));
if compute_stat
    info.nz_blocks = nz_blocks;    % number of zero blocks
    info.rank_block = rank_block;  % rank of each block
    info.norm_block = norm_block;  % nuclear norm of each block
end
info.time = toc(t1);

end


