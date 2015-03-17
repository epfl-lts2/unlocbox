function [sol, info] = proj_nuclearnorm(x, ~, param)
%PROJ_NUCLEARNORM Projection on the nuclear norm ball
%   Usage:  sol=proj_nuclearnorm(x);
%           sol=proj_nuclearnorm(x, gamma, param);
%           [sol,info]=proj_nuclearnorm(...);
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         info  : Structure summarizing informations at convergence
%
%   `proj_nuclearnorm(x, gamma, param)` solves:
%
%   .. sol = argmin_{z} 0.5*||x - z||_2^2  s. t. ||z||_* < epsilon
%
%   .. math::  sol = arg\min_{z} \frac{1}{2} \|x - z\|_2^2 \text{ s. t. } ||z||_* < \epsilon
%   
%   param is a Matlab structure containing the following fields:
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%   * *param.epsilon* : Radius of the nuclear ball (default = 1).
%
%   * *param.svds* : 0 uses svd, 1 uses svds. (default: 1 for sparse
%     matrices, 0 for full matrices)
%
%   * *param.max_rank* : upper bound of rank expected after thresholding.
%     If actual rank is greater, SVDS has to restart with bigger bound.
%     (default: the maximum between 20 and sqrt(n))
%
%   * *param.tol* : tolerance for svds. Bigger tolerance yelds faster
%     results. (default: 1e-5);
%
%   * *param.single* : single precision (1) or not (0)? (default: single
%     only if input is single precision);
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
%   See also:  prox_l1 proj_b1 prox_nuclearnorm

%% TODO: Fix for single precision input! Can't use single and svds at the 
%% same time! Right now "svds" flag overrides "single" flag!


% Authors: Vassilis Kalofolias, Nathanael Perraudin
% Date: Feb 2015 EPFL
%

% Start the time counter
t1 = tic;

if nargin < 3, param = struct; end

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
if param.svds == 1
    if ~isfield(param, 'max_rank'),
        param.max_rank = max(20, sqrt(min(size(x))));
    end
end
if ~isfield(param, 'tol'), param.tol = 1e-5; end
if ~isfield(param, 'single'), param.single = isa(x, 'single'); end
if ~isfield(param, 'epsilon'), param.epsilon = 1; end


svds_opts.tol = param.tol;

if param.svds
    % use upper bound for rank! don't compute everything...
    [U, S, V] = svds(double(x), param.max_rank, 'L', svds_opts);
    %[U, Sigma, V] = svds(x, min(size(x)));
    while sum(S) < param.epsilon
        param.max_rank = 2 * param.max_rank;
        [U, S, V] = svds(double(x), param.max_rank, 'L', svds_opts);
    end
else
    try
        if param.single
            [U, S, V] = svd(single(full(x)), 'econ');       % good for small, dense matrices!!
        else
            [U, S, V] = svd(full(x), 'econ');       % good for small, dense matrices!!
        end
    catch %err
        % This can save you sometimes
        fprintf('SVD failed!! Trying with svds...\n');
        [U, S, V] = svds(double(x), min(size(x)));
    end
end

% Shrink:
sigma = diag(S);            % column vector
paramb1.verbose = 0;
paramb1.epsilon = param.epsilon;
sigma = proj_b1(sigma,0,paramb1);

r = sum(sigma > 0);             % rank of solution

% % This is old code optimized below:
% U = U(:,1:r); V = V(:,1:r); Sigma = diag(sigma(1:r));
% sol = U * Sigma * V';

% Reconstruct X with new singular values
sigma = sigma(1:r);
nuclearNorm = sum(sigma);
% sol = Ur Sr Vr', where Ur = U(:, 1:r), Sr = diag(sigma), Vr = Vr(:, 1:r)
sol = U(:, 1:r) * bsxfun(@times, sigma, V(:, 1:r).');
if param.single
    sol = single(sol);
end

if param.verbose >= 1
    fprintf('  proj nuclear norm: rank= %i, |x|_* = %e \n', r, nuclearNorm);
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


