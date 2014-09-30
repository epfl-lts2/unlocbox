function [sol, info] = prox_tv1d(b, gamma, param)
%PROX_TV1D Total variation proximal operator
%   Usage:  sol=prox_tv1d(x, gamma)
%           sol=prox_tv1d(x, gamma,param)
%           [sol, info]=prox_tv1d(...)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of optional parameters.
%   Output parameters
%         sol   : Solution.
%         info  : Structure summarizing informations at convergence
%
%   This function compute the 1 dimentional TV proximal operator evaluated
%   in b. If b is a matrix, this function will evaluate the TV proximal
%   operator on each row of the matrix. For 2 dimention TV proximal
%   operator the function prox_tv can be used.
%
%   `prox_tv(y, gamma, param)` solves:
%
%   .. sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * ||x||_TV
%
%   .. math::  sol = \min_{z} \frac{1}{2} \|x - z\|_2^2 + \gamma  \|x\|_{TV}
%
%   param is a Matlab structure containing the following fields:
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
%   * *param.useGPU* : Use GPU to compute the TV prox operator. Please prior 
%     call init_gpu and free_gpu to launch and release the GPU library (default: 0).
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
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
%   See also:  prox_tv prox_l1 prox_tv3d 
%
%
%   References: beck2009fastTV


% Author: Nathanael Perraudin, Gilles Puy, Eyal Hirsch
% Date: Jan 2013
%

% Start the time counter
t1 = tic;

% for the GPU
global GLOBAL_useGPU; 

% Optional input arguments

if nargin<3, param=struct; end

if ~isfield(param, 'tol'), param.tol = 10e-4; end
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'useGPU'), param.useGPU = GLOBAL_useGPU; end

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

if param.useGPU
    %gpuDevice(1);
    gamma=gpuArray(gamma);
    if isa(b,'gpuArray')
        allGPU=1;
    else
        b=gpuArray(b);
        allGPU=0;
    end
    % Initializations
    r = gradient_op1d(b*0);
    pold = r; 
    told = gpuArray(1); prev_obj = gpuArray(0); 
    verbose=gpuArray(param.verbose);
    tol=gpuArray(param.tol);
else
    % Initializations
    r = gradient_op1d(b*0);
    pold = r;
    told = 1; prev_obj = 0;
    verbose=param.verbose;
    tol=param.tol;
end

% Main iterations
if verbose > 1
    fprintf('  Proximal TV operator:\n');
end


    
    
for iter = 1:param.maxit

    % Current solution
    sol = b - gamma*div_op1d(r);

    % Objective function value
    obj = .5*norm(b(:)-sol(:), 2)^2 + gamma * sum(norm_tv1d(sol));
    rel_obj = abs(obj-prev_obj)/obj;
    prev_obj = obj;

    % Stopping criterion
    if verbose>1
        fprintf('   Iter %i, obj = %e, rel_obj = %e\n', ...
            iter, obj, rel_obj);
    end
    if rel_obj < tol
        crit = 'TOL_EPS'; break;
    end

    % Udpate divergence vectors and project
    dx = gradient_op1d(sol);
    r = r - 1/(4*gamma) * dx; 
    weights = max(1, abs(r));
    p = r./weights; 

    % FISTA update
    t = (1+sqrt(4*told.^2))/2;
    r = p + (told-1)/t * (p - pold); pold = p;
    told = t;

end


% Log after the minimization
if ~exist('crit', 'var'), crit = 'MAX_IT'; end

if verbose >= 1
    if param.useGPU
        fprintf(['  GPU Prox_TV1D: obj = %e, rel_obj = %e,' ...
            ' %s, iter = %i\n'], obj, rel_obj, crit, iter);
    else
        fprintf(['  Prox_TV1D: obj = %e, rel_obj = %e,' ...
            ' %s, iter = %i\n'], obj, rel_obj, crit, iter);
    end
end



if param.useGPU
    if ~allGPU
        sol=gather(sol);
    end
    info.iter=gather(iter);
    info.final_eval=gather(obj);
else
    info.iter=iter;
    info.final_eval=obj;
end

info.algo=mfilename;
info.crit=crit;
info.time=toc(t1);

end
