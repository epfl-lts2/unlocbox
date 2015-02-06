function [sol, info, objective] = solve_tvdn(y, epsilon, A, At, param)
%SOLVE_TVDN Solve TVDN problem
%   Usage: sol = solve_tvdn(y, epsilon, A, At, param)
%          sol = solve_tvdn(y, epsilon, A, At)
%          [sol,info,objective] = solve_tvdn(...)
%
%   Input parameters:
%         y     : Measurements
%         epsilon: Radius of the L2 ball
%         A     : Operator
%         At    : Adjoint of A
%         param : Optional parameter
%   Output parameters:
%         sol   : Solution
%         info  : Structure summarizing informations at convergence
%         objective: vector (evaluation of the objectiv function each iteration)
%
%   `sol = solve_tvdn(Y, epsilon, A, At, PARAM)` solves:
%
%   .. sol arg min ||x||_TV  s.t.  ||y-A x||_2 < epsilon
%
%   .. math:: arg \min_x \|x\|_{TV}   s.t.  \|y-A x\|_2 < \epsilon
%
%   Y contains the measurements. A is the forward measurement operator and
%   At the associated adjoint operator. PARAM a Matlab structure containing
%   the following fields:
%
%   General parameters:
% 
%   * *param.verbose* : 0 no log, 1 print main steps, 2 print all steps.
%
%   * *param.maxit* : max. nb. of iterations (default: 200).
%
%   * *param.useGPU* : Use GPU to compute the TV prox operator. Please prior 
%     call init_gpu and free_gpu to launch and release the GPU library (default: 0).
%
%   * *param.tol* : is stop criterion for the loop. The algorithm stops if
%
%     ..  (  n(t) - n(t-1) )  / n(t) < tol,
%      
%     .. math:: \frac{  n(t) - n(t-1) }{ n(t)} < tol,
%
%     where  $n(t) = ||(x)||_{TV}$ is the objective function at iteration *t*
%     by default, `tol=10e-4`.
%
%   * *param.gamma* : control the converge speed (default: 1e-1).
% 
% 
%   Projection onto the L2-ball :
%
%   * *param.tight_b2* : 1 if A is a tight frame or 0 if not (default = 1)
% 
%   * *param.nu_b2* : bound on the norm of the operator A, i.e.
%
%     .. ` ||A x||^2 <= nu * ||x||^2 
%
%     .. math::  \|A x\|^2 \leq \nu  \|x\|^2 
%
%   * *param.tol_b2* : tolerance for the projection onto the L2 ball (default: 1e-3):
%
%     .. epsilon/(1-tol) <= ||y - A z||_2 <= epsilon/(1+tol)
%
%     .. math:: \frac{\epsilon}{1-tol} \leq \|y - A z\|_2 \leq \frac{\epsilon}{1+tol}
%    
%   * *param.maxit_b2* : max. nb. of iterations for the projection onto the L2
%     ball (default 200).
% 
% 
%   Proximal TV operator:
%
%   * *param.maxit_tv* : Used as stopping criterion for the proximal TV
%     operator. Maximum number of iterations.
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
%   * *info.final_eval* : Final evaluation of the objectivs functions
%
%   * *info.crit* : Stopping critterion used 
%
%   * *info.rel_norm* : Relative norm at convergence 
%
%   * *info.residue* : Final residue 
%
%
%   The problem is solved thanks to a Douglas-Rachford splitting
%   algorithm.
%
%   Demos: demo_tvdn
%
%   References: combettes2007douglas

% Author: Gilles Puy, Nathanael Perraudin
% Date: Nov. 1, 2012

% Start the time counter
t1 = tic;

% Optional input arguments
if nargin<5, param=struct; end


% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'tol'), param.tol = 1e-4; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'gamma'), param.gamma = 1e-2; end
if ~isfield(param, 'useGPU'), param.useGPU = 0; end

% Input arguments for projection onto the L2 ball
param_b2.A = A; param_b2.At = At;
param_b2.y = y; param_b2.epsilon = epsilon;
param_b2.verbose = param.verbose;
if isfield(param, 'nu_b2'), param_b2.nu = param.nu_b2; end
if isfield(param, 'tol_b2'), param_b2.tol = param.tol_b2; end
if isfield(param, 'tight_b2'), param_b2.tight = param.tight_b2; end
if isfield(param, 'maxit_b2')
    param_b2.maxit = param.maxit_b2;
end

% Input arguments for prox TV
param_tv.verbose = param.verbose; param_tv.tol = param.tol;
param_tv.useGPU = param.useGPU;
if isfield(param, 'maxit_tv')
    param_tv.maxit= param.maxit_tv;
end

% Initialization
xhat = At(y); 
[~,~,prev_norm,iter,objective,~] = convergence_test(norm_tv(y));
% Main loop
while 1
    
    %
    if param.verbose>=1
        fprintf('Iteration %i:\n', iter);
    end
    
    % Projection onto the L2-ball
    [sol, param_b2.u] = proj_b2(xhat, NaN, param_b2);
    
    % Global stopping criterion
    curr_norm = norm_tv(sol);
    [stop,rel_norm,prev_norm,iter,objective,crit] =...
        convergence_test(curr_norm,prev_norm,iter,objective,param);
    [sol,param] = post_process(sol, iter, curr_norm, prev_norm, objective, param);
    if stop
        break;
    end
    if param.verbose >= 1
        fprintf('  ||x||_TV = %e, rel_norm = %e\n', ...
            curr_norm, rel_norm);
    end

    
    % Proximal L1 operator
    xhat = 2*sol - xhat;
    temp = prox_tv(xhat, param.gamma, param_tv);
    xhat = temp + sol - xhat;
    
    
end

% Residual
temp = A(sol);
res=norm(y(:)-temp(:));
% Log
if param.verbose>=1
    % L1 norm
    fprintf('\n Solution found:\n');
    fprintf(' Final TV norm: %e\n', curr_norm);
    
    % Residual
    
    fprintf(' epsilon = %e, ||y-Ax||_2=%e\n', epsilon, ...
        res);
    
    % Stopping criterion
    fprintf(' %i iterations\n', iter);
    fprintf(' Stopping criterion: %s \n\n', crit);
    
end


info.algo=mfilename;
info.iter=iter;
info.final_eval=curr_norm;
info.crit=crit;
info.time=toc(t1);
info.rel_norm=rel_norm;
info.residue=res;

end
