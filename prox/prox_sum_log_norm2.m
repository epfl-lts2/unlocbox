function [sol,info] = prox_sum_log_norm2(x, alpha, beta, gamma, param)
%PROX_SUM_LOG_NORM2 Proximal operator of log-barrier  - sum(log(x))
%   Usage:  sol = prox_sum_log_norm2(x, alpha, beta, gamma)
%           sol = prox_sum_log_norm2(x, alpha, beta, gamma, param)
%           [sol, info] = prox_sum_log_norm2(x, alpha, beta, gamma, param)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         alpha : multiplier of -log
%         beta  : multiplier of norm-2
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         info : Structure summarizing informations at convergence
%
%   `prox_l1(x, gamma, param)` solves:
%
%   .. sol = argmin_{z} 0.5*||x - z||_2^2 - gamma * (alpha*sum(log(z)) + beta/2 ||z||^2)
%
%   .. math::  sol = \arg\min_{z} \frac{1}{2} \|x - z\|_2^2 - \gamma (\alpha \sum_i(log(z_i)) + \beta/2 \|z\|_2^2)
%
%   param is a Matlab structure containing the following fields:
%
%   * *param.verbose* : 0 no log, (1) print -sum(log(z)), 2 additionally
%       report negative inputs.
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
%   * *info.crit* : Stopping critterion used 
%
%
%   See also:  prox_sum_log
%

% Author: Vassilis Kalofolias
% Date: September 2015
% Testing: TODO test_prox_sum_log_norm2

% Start the time counter
t1 = tic;

if nargin < 3, param = struct; end

if ~isfield(param, 'verbose'), param.verbose = 1; end


% test the parameters
% stop_error = test_gamma(gamma);
if gamma < 0
    error('Gamma can not be negative');
elseif gamma == 0
    stop_error = 1;
else
    stop_error = 0;
end

if stop_error
    sol = x;
    info.algo = mfilename;
    info.iter = 0;
    info.final_eval = 0;
    info.crit = '--';
    info.time = toc(t1);
    return
end



sol = (x + sqrt(x.^2 + 4*alpha*gamma*(1+beta*gamma))) /(2*(1+beta*gamma));
info.algo = mfilename;
info.iter = 0;
info.final_eval = gamma * (-alpha*sum(log(x(:))) + beta*norm(x(:))^2);
info.crit = '--';
info.time = toc(t1);
    
% Log after the prox
if param.verbose >= 1
    fprintf('  prox_sum_log: - sum(log(x)) = %e', info.final_eval / gamma);
    if param.verbose > 1
        n_neg = nnz(x(:) <= 0);
        if n_neg > 0
            fprintf(' (%d negative elements, log not defined, check stability)', n_neg);
        end
    end
    fprintf('\n');
end


end


