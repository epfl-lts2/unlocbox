function [sol, info] = prox_nuclearnorm(x, gamma, param)
%PROX_NUCLEARNORM Proximal operator with the nuclear norm
%   Usage:  sol=prox_nuclearnorm(x, gamma)
%           sol=prox_nuclearnorm(x, gamma, param)
%           [sol,info]=prox_nuclearnorm(...)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         info  : Structure summarizing informations at convergence
%
%   `prox_NuclearNorm(x, gamma, param)` solves:
%
%   .. sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * ||z||_*
%
%   .. math::  sol = arg\min_{z} \frac{1}{2} \|x - z\|_2^2 + \gamma  ||z||_*
%   
%   param is a Matlab structure containing the following fields:
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%   * *param.svds* : 0 uses svd, 1 uses svds. (default: 1 for sparse
%     matrices, 0 for full matrices)
%
%   * *param.max_rank* : TODO!! will be availlable soon for large sparse
%     matrices and svds!! 
%
%   info is a Matlab structure containing the following fields:
%
%   * *info.algo* : Algorithm used
%
%   * *info.iter* : Number of iteration
%
%   * *info.time* : Time of exectution of the function in sec.
%
%   * *info.final_eval* : Final evaluation of the function
%
%   * *info.crit* : Stopping criterion used 
%
%   * *info.rank* : Rank of the final solution (-1 means the rank was not
%     computed) 
%
%
%   See also:  prox_l1 proj_b1 prox_tv


% Authors: Nathanael Perraudin, Vassilis Kalofolias
% Date: June 2012 EPFL
%

% Start the time counter
t1 = tic;

if nargin<3, param=struct; end

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'svds'), 
    if issparse(x)
        % use svds for sparse matrices bigger than 1000 x 1000
        param.svds = 1;
    else
        % otherwise use svd (econ mode)
        param.svds = 0;     
    end
end


% Test of gamma
if test_gamma(gamma)
    sol = x;
    % set the information structure returned
    iter = 0;
    crit = '--';
    
    info.algo = mfilename;
    info.iter = iter;
    info.final_eval = 0;
    info.crit = crit;
    info.time = toc(t1);
    info.rank = -1; %rank(full(sol));
    return
end

% Useful functions
soft = @(z, T) sign(z).*max(abs(z) - T, 0);

if param.svds
    %TODO: use upper bound for rank! don't compute everything...
    [U, Sigma, V] = svds(x, min(size(x)));
else
    try
        [U, Sigma, V] = svd(full(x), 'econ');       % good for small, dense matrices!!
    catch %err
        % This can save you sometimes
        fprintf('SVD failed!! Trying with svds...\n');
        [U, Sigma, V] = svds(x, min(size(x)));
    end
end

% Shrink:
sigma = diag(Sigma);            % column vector
sigma = soft(sigma, gamma);     % modified singular values!
r = sum(sigma > 0);             % rank of solution

% % This is old code optimized below:
% U = U(:,1:r); V = V(:,1:r); Sigma = diag(sigma(1:r));
% sol = U * Sigma * V';

% Reconstruct X with new singular values
sigma = sigma(1:r);
nuclearNorm = sum(sigma);
% sol = Ur Sr Vr', where Ur = U(:, 1:r), Sr = diag(sigma), Vr = Vr(:, 1:r)
sol = U(:, 1:r) * bsxfun(@times, sigma, V(:, 1:r).');

if param.verbose >= 1
    fprintf('  prox nuclear norm: rank= %i, |x|_* = %e \n', r, nuclearNorm);
end

% set the information structure returned
iter = 0;
crit = '--';
info.algo = mfilename;
info.iter = iter;
info.final_eval = nuclearNorm;
info.crit = crit;
info.time = toc(t1);
info.rank = r;

end


