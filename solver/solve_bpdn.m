function [sol, info] = solve_bpdn(y, epsilon, A, At, Psi, Psit, param)
%SOLVE_BPDN Solve BPDN (basis pursuit denoising) problem
%   Usage: sol = solve_bpdn(y, epsilon, A, At, Psi, Psit, param)
%          sol = solve_bpdn(y, epsilon, A, At, Psi, Psit)
%          [sol, info] = solve_bpdn(...)
%
%   Input parameters:
%         y     : Measurements
%         epsilon: Radius of the L2 ball
%         A     : Operator
%         At    : Adjoint of A
%         Psi   : Operator
%         Psit  : Adjoint of Psi
%         param : Optional parameter
%   Output parameters:
%         sol   : Solution
%         info  : Structure summarizing informations at convergence
%
%   `sol = solve_BPDN(y, A, At, Psi, Psit, param)` solves:
%
%   .. sol arg min ||Psi x||_1   s.t.  ||y-A x||_2 < epsilon
%
%   .. math:: arg \min_x \| \Psi x\|_1   s.t.  \|y-A x\|_2 < \epsilon
%
%   Y contains the measurements. A is the forward measurement operator and
%   At the associated adjoint operator. Psit is a sparfying transform and Psi
%   its adjoint. PARAM a Matlab structure containing the following fields:
%
%   General parameters:
% 
%   * *param.verbose* : 0 no log, 1 print main steps, 2 print all steps.
%
%   * *param.maxit* : max. nb. of iterations (default: 200).
%
%   * *param.tol* : is stop criterion for the loop. The algorithm stops if
%
%     ..  (  n(t) - n(t-1) )  / n(t) < tol,
%      
%     .. math:: \frac{  n(t) - n(t-1) }{ n(t)} < tol,
%
%     where  $n(t) = ||\Psi(x)||$ is the objective function at iteration *t*
%     by default, `tol=10e-4`.
%
%   * *param.gamma* : control the converge speed (default: 1).
% 
% 
%   Projection onto the L2-ball :
%
%   * *param.tight_b2* : 1 if A is a tight frame or 0 if not (default = 1)
% 
%   * *nu_b2* : bound on the norm of the operator A, i.e.
%
%     .. ` ||A x||^2 <= nu * ||x||^2 
%
%     .. math::  \|A x\|^2 \leq \nu  \|x\|^2 
%
%   * *tol_b2* : tolerance for the projection onto the L2 ball (default: 1e-3):
%
%   .. epsilon/(1-tol) <= ||y - A z||_2 <= epsilon/(1+tol)
%
%   .. math:: \frac{\epsilon}{1-tol} \leq \|y - A z\|_2 \leq \frac{\epsilon}{1+tol}
%    
%   * *maxit_b2* : max. nb. of iterations for the projection onto the L2
%     ball (default 200).
% 
% 
%   Proximal L1 operator:
%
%   * *tol_l1* : Used as stopping criterion for the proximal L1
%     operator. Min. relative change of the objective value between two
%     successive estimates.
%
%   * *maxit_l1* : Used as stopping criterion for the proximal L1
%     operator. Maximum number of iterations.
% 
%   * *param.nu_l1* : bound on the norm^2 of the operator Psi, i.e.
%
%     .. ` ||Psi x||^2 <= nu * ||x||^2 
%
%     .. math::  \|\Psi x\|^2 \leq \nu  \|x\|^2 
% 
%   * *param.tight_l1* : 1 if Psit is a tight frame or 0 if not (default = 1)
% 
%   * *param.weights* : weights (default = 1) for a weighted L1-norm defined
%     as:
%
%     .. sum_i{weights_i.*abs(x_i)}
%
%     .. math::  \sum_i w_i |x_i|     
%
%   The problem is solved thanks to a Douglas-Rachford splitting
%   algorithm.
%
%   Demos: demo_weighted_l1
%
%   References: combettes2007douglas

% Author: Gilles Puy, Nathanael Perraudin
% Date: Nov. 1, 2012
%


% Optional input arguments
if nargin<7, param=struct; end
if nargin < 5; Psi = @(x) x; end
if nargin < 6; Psit = Psi; end

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'tol'), param.tol = 1e-4; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'pos_l1'), param.pos_l1 = 0; end

% Input arguments for projection onto the L2 ball
param_b2.A = A; param_b2.At = At;
param_b2.y = y; param_b2.epsilon = epsilon;
param_b2.verbose = param.verbose-1;
if isfield(param, 'nu_b2'), param_b2.nu = param.nu_b2; end
if isfield(param, 'tol_b2'), param_b2.tol = param.tol_b2; end
if isfield(param, 'tight_b2'), param_b2.tight = param.tight_b2; end
if isfield(param, 'maxit_b2')
    param_b2.maxit = param.maxit_b2;
end


% Input arguments for prox L1
param_l1.A = Psi; param_l1.At = Psit; param_l1.pos = param.pos_l1;
param_l1.verbose = param.verbose-1; param_l1.tol = param.tol;
if isfield(param, 'nu_l1')
    param_l1.nu = param.nu_l1;
end
if isfield(param, 'tight_l1')
    param_l1.tight = param.tight_l1;
end
if isfield(param, 'maxit_l1')
    param_l1.maxit = param.maxit_l1;
end
if isfield(param, 'tol_l1')
    param_l1.tol = param.tol_l1;
end
if isfield(param, 'weights')
    param_l1.weights = param.weights;
    f1.eval = @(x) norm(param_l1.weights(:).*x(:),1);
else
    param_l1.weights = 1;
    f1.eval = @(x) norm(x(:),1);
end

f1.prox = @(x,T) prox_l1(x,T,param_l1);


f2.prox = @(x,T) proj_b2(x,T,param_b2);
f2.eval = @(x) eps;

[sol, info] = douglas_rachford(At(y),f2, f1, param);

end
