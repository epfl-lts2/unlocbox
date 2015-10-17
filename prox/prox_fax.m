function [sol, info]  = prox_fax( z,gamma,param )
%PROX_FAX Proximal operator of the adjoint function of f
%   Usage:   sol=prox_adjoint(x, gamma, param);
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter (ususally it should be 1)
%         param : Parameter (Please see: param.f)
%   Output parameters:
%         sol   : Solution.
%         infos : Structure summarizing informations at convergence
%
%   ` prox_fax( x,gamma,param )` solves:
%
%   .. sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * f(Ax)
%
%   .. math::  sol = \min_{z} \frac{1}{2} \|x - z\|_2^2 + \gamma  f(Ax)
%
%   This method allows to compute the proximal operator of f(Ax) when only
%   the proximal operator of A can be computed. This function use an ADMM
%   splitting.
%
%   *param* is a non optional structure of parameter containing 2
%   mendatory parameter:
%
%   * *param.A* : Forward operator 
%
%   * *param.At* : Adjoint operator
%
%   * *param.f* : is a structure representing a convex function.
%     Inside the structure, there have to be the prox of the function that
%     can be called by *f1.prox* and the function itself that can be called
%     by *f1.eval*.  
%
%   As an option, you may specify
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
%   * *param.L2_maxit* : max. nb. of iterations for the l2 proximal
%     operator (default: 30). 
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)


% Author: Nathanael Perraudin
% Date: 31 May 2013
% Testing: test_prox_l1

if nargin<3
   error('PROX_FAX: you need to provide param for this function')
end

if ~isfield(param,'A') || ~isfield(param,'At')
    error('PROX_FAX: the field A and At of param are mendatory')
end

if ~isfield(param,'f') 
    error('PROX_FAX: the field f of param is mendatory')
end

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'tight'), param.tight = 1; end
if ~isfield(param, 'nu'), param.nu = 1; end
if ~isfield(param, 'tol'), param.tol = 1e-3; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'maxit_l2'), param.maxit_l2 = 30; end


f1.eval = @(x) norm(reshape(x,[],1))^2;
paraml2.A = param.A;
paraml2.At = param.At;
paraml2.maxit = param.maxit_l2;
paraml2.verbose = param.verbose - 1;
paraml2.tol = param.tol;
paraml2.nu = param.nu;
paraml2.tight = param.tight;
f1.proxL = @(x,T) reverse_prox(x,0.5*T,paraml2,z);

f2.eval = @(x) gamma*param.f.eval(x);
f2.prox = @(x,T) param.f.prox(x,T*gamma);
f2.L = param.A;

paramsolver.maxit = param.maxit;
paramsolver.verbose = param.verbose;
paramsolver.tol = param.tol;

[sol, info] = admm(z,f1, f2, paramsolver);


end

function sol = reverse_prox(x,T,paraml2,z)
    paraml2.y = x;
    sol = prox_l2(z,T,paraml2);
end
