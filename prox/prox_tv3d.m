function [sol,info] = prox_tv3d(x, gamma, param)
%PROX_TV3D Total variation proximal operator
%   Usage:  sol=prox_tv3d(x, gamma)
%           sol=prox_tv3d(x, gamma,param)
%           [sol, info]=prox_tv3d(...)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         info  : Structure summarizing informations at convergence
%
%   This function compute the 3 dimentional TV proximal operator evaluated
%   in b. If b is 4 dimentional, this function will evaluate the TV
%   proximal operator on each cube. For 2 dimention TV proximal of cubes
%   operator the function prox_tv can be used.
%
%   `prox_tv3d(y, gamma, param)` solves:
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
%   * *param.parrallel* : Parallelisation level. 0 means no
%     parallelization, 1 means all cubes (fourth dimension changing) at the
%     same time.
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%   * *param.useGPU* : Use GPU to compute the TV prox operator. Please prior 
%     call init_gpu and free_gpu to launch and release the GPU library (default: 0).
%
%   * *param.weights* : weights for each dimention (default $[1, 1, 1]$)
%
%   infos is a Matlab structure containing the following fields:
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
%   See also:  prox_l1 prox_tv
%
%   References: beck2009fastTV


% Author: Nathanael Perraudin, William Guicquero
% Date: October 15, 2010
%

% Start the time counter
t1 = tic;

% for the GPU
global GLOBAL_useGPU; 

% Optional input arguments
if nargin<3, param=struct; end

if ~isfield(param, 'tol'), param.tol = 1e-4; end
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'useGPU')
    param.useGPU = (GLOBAL_useGPU) || (isa(x,'gpuArray')); 
end
if ~isfield(param, 'weights'), param.weights = [1,1,1]; end


if ~isfield(param, 'parallel')
    if size(x,4)==1
        param.parallel = 1;
    else
        param.parallel = 0; 
    end
end



if param.parallel == 0
   % call prox 3d for each cube
   param.parallel = 1;
   sol = zeros(size(x));
   info.iter = 0;
   info.time = 0;
   info.algo=mfilename;
   info.final_eval = 0;
   info.crit = 'TOL_EPS';  % return this only if ALL subproblems finish with this criterion.
   param.verbose = param.verbose-1; % Handle verbosity
   
   for ii = 1:size(x, 4)
       [sol(:, :, :, ii), infos_ii] = prox_tv3d(x(:,:,:,ii), gamma, param);
       info.iter = info.iter + infos_ii.iter;
       info.time = info.time + infos_ii.time;
       info.final_eval = info.final_eval + infos_ii.final_eval;

       if strcmpi(infos_ii.crit, 'MAX_IT');
           info.crit = 'MAX_IT';   % if ANY subproblem reaches maximum iterations, return this as criterion!
       end
   end
   
   return
   
end

% If once parfor is working generally on MATLAB
% if strcmpi(param.parallel, 'parfor')
%    % recall prox 3d for each cube
%    param.parallel = 'full';
%    sol = zeros(size(x));
%    
%    parfor ii = 1:size(x,4)
%        sol(:,:,:,ii) = prox_tv3d(x(:,:,:,ii), gamma, param);
%    end
%    
%    return
%    
% end


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


wx = param.weights(1);
wy = param.weights(2);
wz = param.weights(3);
mt = max(param.weights);

% Initializations
if param.useGPU
    %gpuDevice(1);
    gamma=gpuArray(gamma);
    if isa(x,'gpuArray')
        allGPU=1;
    else
        x=gpuArray(x);
        allGPU=0;
    end
    % Initializations
    [r, s, k] = gradient_op3d(x*0);
    pold = r; qold = s; kold = k;
    told = gpuArray(1); prev_obj = gpuArray(0); 
    verbose=gpuArray(param.verbose);
    tol=gpuArray(param.tol);
else
    [r, s, k] = gradient_op3d(x*0);
    pold = r; qold = s; kold = k;
    told = 1; prev_obj = 0;
    verbose=param.verbose;
    tol=param.tol;
end

% Main iterations
if verbose > 1
    if param.useGPU
        fprintf('  Proximal TV operator using TV:\n');
    else
        fprintf('  Proximal TV operator:\n');
    end
end


for iter = 1:param.maxit
    
    % Current solution
    sol = x - gamma * div_op3d(r, s, k,wx,wy,wz);
    
    % Objective function value
    obj = .5*norm(x(:)-sol(:), 2)^2 + gamma * sum(norm_tv3d(sol,wx,wy,wz));
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
    % TODO: read reference for good explanation... We change lemma 4.2 to
    % be valid for 3D denoising and we should get a bound with 12 instead
    % of 8.
    [dx, dy, dz] = gradient_op3d(sol,wx,wy,wz);
    r = r - 1/(12*gamma*mt^2) * dx;
    s = s - 1/(12*gamma*mt^2) * dy;
    k = k - 1/(12*gamma*mt^2) * dz;
    % Isotropic tv
    weights = max(1, sqrt(abs(r).^2+abs(s).^2+abs(k).^2));
    % anisotropic TV
    %weights = max(1, abs(r)+abs(s)+abs(k));
    p = r./weights;
    q = s./weights;
    o = k./weights;
    
    
    % FISTA update
    t = (1+sqrt(4*told^2))/2;
    r = p + (told-1)/t * (p - pold); pold = p;
    s = q + (told-1)/t * (q - qold); qold = q;
    k = o + (told-1)/t * (o - kold); kold = o;
    told = t;
    
end

% Log after the minimization
if ~exist('crit', 'var'), crit = 'MAX_IT'; end



if verbose >= 1
    if param.useGPU
        fprintf(['  GPU Prox_TV 3D: obj = %e, rel_obj = %e,' ...
            ' %s, iter = %i\n'], obj, rel_obj, crit, iter);
    else
        fprintf(['  Prox_TV 3D: obj = %e, rel_obj = %e,' ...
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
info.iter=iter;
info.final_eval=obj;
info.crit=crit;
info.time=toc(t1);

end
