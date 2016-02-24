function [sol,info] = prox_l1_admm(z, gamma, param)
%PROX_L1_ADMM Proximal operator with L1 norm
%   Usage:  sol=prox_l1_admm(x, gamma)
%           sol=prox_l1_admm(x, gamma, param)
%           [sol, info]=prox_l1_admm(x, gamma, param)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         info : Structure summarizing informations at convergence
%
%   `prox_l1(x, gamma, param)` solves:
%
%   .. sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * ||A ( z - y )||_1
%
%   .. math::  sol = \min_{z} \frac{1}{2} \|x - z\|_2^2 + \gamma  \|A z\|_1
%
%   param is a Matlab structure containing the following fields:
%
%   * *param.A* : Forward operator (default: Id).
%
%   * *param.At* : Adjoint operator (default: Id).
%
%   * *param.y* : (Under developement).
%
%   * *param.tight* : 1 if A is a tight frame or 0 if not (default = 1)
%
%   * *param.nu* : bound on the norm of the operator A (default: 1), i.e.
%
%     .. ` ||A x||^2 <= nu * ||x||^2 
%
%     .. math::  \|A x\|^2 \leq \nu  \|x\|^2 
%   
%   * *param.tol* : is stop criterion for the loop. The algorithm stops if
%
%     ..  (  n(t) - n(t-1) )  / n(t) < tol,
%      
%     .. math:: \frac{  n(t) - n(t-1) }{ n(t)} < tol,
%
%     where  $n(t) = f(x)+ 0.5 \|x-z\|_2^2$ is the objective function at iteration *t*
%     by default, `tol=10e-4`.
%
%   * *param.maxit* : max. nb. of iterations (default: 200).
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%   * *param.weights* : weights for a weighted L1-norm (default = 1)
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
%   See also:  proj_b1 prox_linf1 prox_l21 prox_tv
%
%   References: fadili2009monotone boyd2011distributed van2008probing

% Author: Nathanael Perraudin
% Date: Nov 2012
%

% Start the time counter
t1 = tic;

% Optional input arguments
if nargin<3, param=struct; end

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'tight'), param.tight = 1; end
if ~isfield(param, 'nu'), param.nu = 1; end
if ~isfield(param, 'tol'), param.tol = 1e-3; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'At'), param.At = @(x) x; end
if ~isfield(param, 'A'), param.A = @(x) x; end
if ~isfield(param, 'y'), param.y = zeros(size(z)); end
if ~isfield(param, 'weights'), param.weights = 1; end
if ~isfield(param, 'maxit_l2'), param.maxit_l2 = 30; end

% test the parameters
if test_gamma(gamma)
    sol = x;
    info.algo=mfilename;
    info.iter=0;
    info.final_eval=0;
    info.crit='--';
    info.time=toc(t1);
    return; 
end

param.weights = test_weights(param.weights);


f1.eval = @(x) norm(reshape(x,[],1))^2;
paraml2.A = param.A;
paraml2.At = param.At;
paraml2.maxit = param.maxit_l2;
paraml2.verbose = param.verbose - 1;
paraml2.tol = param.tol;
paraml2.nu = param.nu;
paraml2.tight = param.tight;
f1.proxL = @(x,T) reverse_prox(x,0.5*T,paraml2,z);

f2.eval = @(x) gamma*norm(x-param.y,1);
y = param.y;
paraml1.verbose = param.verbose -1;
f2.prox = @(x,T) prox_l1(x-y,T*gamma,paraml1)+y;
f2.L = param.A;

paramsolver.maxit = param.maxit;
paramsolver.verbose = param.verbose;
paramsolver.tol = param.tol;

[sol, info] = admm(z,f1, f2, paramsolver);


info.algo=mfilename;




end

function sol = reverse_prox(x,T,paraml2,z)
    paraml2.y = x;
    sol = prox_l2(z,T,paraml2);
end

