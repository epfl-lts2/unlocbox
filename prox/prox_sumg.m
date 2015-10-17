function [sol, info] = prox_sumg(x, gamma , param)
%PROX_sumG Proximal operator of a sum of function
%   Usage:  sol=prox_sumg(x, gamma, param)
%           [sol, info]=prox_sumg(...)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of parameters.
%   Output parameters:
%         sol   : Solution.
%         info  : Structure summarizing informations at convergence
%
%   `prox_sumG(x, gamma , param)` solves:
%
%   .. sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * Sum_{i} w_i G_i(z)     for z,x belong to R^N
%
%   .. math::  sol = \operatorname*{arg\,min}_{z} \frac{1}{2} \|x - z\|_2^2 + \gamma  \sum_{i} w_i G_i(z)      \hspace{1cm} for \hspace{1cm}  z,x\in R^N
%
%   param is a Matlab structure containing the following fields:
%
%   * *param.G* : cellarray of structure with all the prox operator inside and eventually 
%     the norm if no norm is defined, the $L^1$ norm is used the prox: *F{i}.prox* and 
%     norm: *F{i}.eval* are defined in the same way as in the
%     Forward-backward and Douglas-Rachford algorithms
%
%   * *param.weights* : weights of different functions (default = $1/N$,
%     where $N$ is the total number of function) 
%
%   * *param.tol* : is stop criterion for the loop. The algorithm stops if
%
%     ..  (  n(t) - n(t-1) )  / n(t) < tol,
%      
%     .. math:: \frac{  n(t) - n(t-1) }{ n(t)} < tol,
%
%     where  $n(t) = f_1(Lx)+f_2(x)$ is the objective function at iteration *t*
%     by default, `tol=10e-4`.
%
%   * *param.lambda_t*: is the weight of the update term. By default 1.
%     This should be between 0 and 1.
%
%   * *param.maxit* : is the maximum number of iteration. By default, it is 200.
% 
%   * *param.verbose* : 0 no log, 1 print main steps, 2 print all steps.     
%
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
%   See also:  generalized_forward_backward
%
%   Demo: demo_prox_multi_functions
%
%   References: raguet2011generalized


% Author:  Nathanael Perraudin
% Date: Mars 2015
% Testing: test_prox_functions

% Optional input arguments
if nargin<3, param=struct; end


% Test of gamma
if test_gamma(gamma)
    sol = x;
    info.algo=mfilename;
    info.iter=0;
    info.final_eval=0;
    info.crit='--';
    info.time=toc(t1);
    return; 
end

if ~isfield(param, 'G'), param.G = cell(0); end
% Definition of the gradient function
f.grad = @(y) 1/gamma * (y-x);
f.eval = @(y) 0.5/gamma * norm(y-x,'fro')^2;
f.beta = 1/gamma;

[sol, info] = solvep(x,{f,param.G{1:end}},param);

end
