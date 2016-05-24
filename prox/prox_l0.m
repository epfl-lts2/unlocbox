function [sol, info] = prox_l0(x, gamma, param)
%PROX_L0 Proximal operator of the L0 norm
%   Usage:  sol = prox_l0(x)
%           sol = prox_l0(x, gamma)
%           sol = prox_l0(x, gamma, param)
%           [sol, info] = prox_l0(x, gamma, param)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         info  : Structure summarizing information at convergence
%
%   `prox_l0(x, gamma, param)` solves:
%
%   .. sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * || z ||_0
%
%   .. math::  sol = \min_{z} \frac{1}{2} \|x - z\|_2^2 + \gamma  \| z\|_0
%
%   param is a Matlab structure containing the following fields:
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%   * *param.k* : number of non zero elements (if not defined, it uses
%     gamma to dertermine how many coefficients are kept
%
%   info is a Matlab structure containing the following fields:
%
%   * *info.algo* : Algorithm used
%
%   * *info.iter* : Number of iteration
%
%   * *info.time* : Time of execution of the function in sec.
%
%   * *info.final_eval* : Final evaluation of the function
%
%   * *info.crit* : Stopping criterion used 
%
%   See also:  proj_b1 prox_l1


% Author: Nathanael Perraudin
% Date: Nov 2012
% Testing: test_prox_l1

% Start the time counter
t1 = tic;

% Optional input arguments
if nargin<3, param = struct; end

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end

% test the parameters
if test_gamma(gamma)
    sol = x;
    info.algo = mfilename;
    info.iter = 0;
    info.final_eval = 0;
    info.crit = '--';
    info.time = toc(t1);
    return; 
end

if isfield(param, 'k') % k defined
    [~, ind] = sort(x, 'ascend');
    sol = x;
    sol(ind(1:(end-param.k))) = 0;
    norm_l0 = param.k;    
else % k not defined
    hard = @(x,T) x .* (abs(x) > T);
    sol = hard(x, gamma);
    norm_l0 = sum(abs(sol(:)) > 0);
end

crit = 'TOL_EPS';
iter = 1;

% Print log
if param.verbose >= 1
    fprintf(['  prox_L0: ||x||_0 = %e,', ...
        ' %s, iter = %i\n'], norm_l0, crit, iter);
end

info.algo = mfilename;
info.iter = iter;
info.final_eval = norm_l0;
info.crit = crit;
info.time = toc(t1);

end


