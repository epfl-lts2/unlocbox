function [sol,info] = prox_tv4d(x, gamma, param)
%PROX_TV4D Total variation proximal operator
%   Usage:  sol=prox_tv4d(x, gamma)
%           sol=prox_tv4d(x, gamma,param)
%           [sol, info]=prox_tv4d(...)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         info  : Structure summarizing informations at convergence
%
%   This function compute the 4 dimentional TV proximal operator evaluated
%   in b. If b is 5 dimentional, this function will evaluate the TV
%   proximal operator on each 4 dimentional cube. 
%
%   `prox_tv4d(y, gamma, param)` solves:
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
%     parallelization, 1 means working on all the data at once
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%   * *param.weights* : weights for each dimention (default $[1, 1, 1, 1]$)
%
%   * *param.useGPU* : Use GPU to compute the TV prox operator. Please prior 
%     call init_gpu and free_gpu to launch and release the GPU library (default: 0).
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
% Date: April 25, 2014
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
if ~isfield(param, 'useGPU'), param.useGPU = GLOBAL_useGPU; end
if ~isfield(param, 'weights'), param.weights = [1,1,1,1]; end


if ~isfield(param, 'parallel')
    if size(x,5)==1
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
   
   for ii = 1:size(x, 5)
       [sol(:, :, :,:, ii), infos_ii] = prox_tv4d(x(:,:,:,:,ii), gamma, param);
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
%        sol(:,:,:,:,ii) = prox_tv3d(x(:,:,:,:,ii), gamma, param);
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
wt = param.weights(4);
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
    [r, s, k, u] = gradient_op4d(x*0);
    pold = r; qold = s; kold = k; uold = u;
    told = gpuArray(1); prev_obj = gpuArray(0); 
    verbose=gpuArray(param.verbose);
    tol=gpuArray(param.tol);
else
    [r, s, k, u] = gradient_op4d(x*0);
    pold = r; qold = s; kold = k; uold = u;
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
    sol = x - gamma * div_op4d(r, s, k, u, wx, wy, wz, wt);
    
    % Objective function value
    obj = .5*norm(x(:)-sol(:), 2)^2 + ...
        gamma * sum(norm_tv4d(sol, wx, wy, wz, wt));
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
    % be valid for 3D denoising and we should get a bound with 16 instead
    % of 8.
    [dx, dy, dz, dt] = gradient_op4d(sol, wx, wy, wz, wt);
    r = r - 1/(16*gamma*mt^2) * dx;
    s = s - 1/(16*gamma*mt^2) * dy;
    k = k - 1/(16*gamma*mt^2) * dz;
    u = u - 1/(16*gamma*mt^2) * dt;

    % Isotropic tv
    weights = max(1, sqrt(abs(r).^2+abs(s).^2+abs(k).^2+abs(u).^2));
    % anisotropic TV
    %weights = max(1, abs(r)+abs(s)+abs(k));
    p = r./weights;
    q = s./weights;
    o = k./weights;
    m = u./weights;

    
    
    % FISTA update
    t = (1+sqrt(4*told^2))/2;
    r = p + (told-1)/t * (p - pold); pold = p;
    s = q + (told-1)/t * (q - qold); qold = q;
    k = o + (told-1)/t * (o - kold); kold = o;
    u = m + (told-1)/t * (m - uold); uold = m;

    told = t;
    
end

% Log after the minimization
if ~exist('crit', 'var'), crit = 'MAX_IT'; end



if verbose >= 1
    if param.useGPU
        fprintf(['  GPU Prox_TV 4D: obj = %e, rel_obj = %e,' ...
            ' %s, iter = %i\n'], obj, rel_obj, crit, iter);
    else
        fprintf(['  Prox_TV 4D: obj = %e, rel_obj = %e,' ...
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
