function [sol, info] = prox_tv1d(x, gamma, param)
%PROX_TV1D Total variation proximal operator
%   Usage:  sol=prox_tv1d(x, gamma)
%           sol=prox_tv1d(x, gamma, param)
%           [sol, info]=prox_tv1d(...)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of optional parameters.
%   Output parameters
%         sol   : Solution.
%         info  : Structure summarizing information at convergence
%
%   This function computes the 1 dimensional TV proximal operator evaluated
%   in b. If b is a matrix, this function will evaluate the TV proximal
%   operator on each row of the matrix. For 2D, TV proximal
%   operator `prox_tv` can be used.
%
%   `prox_tv(y, gamma, param)` solves:
%
%   .. sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * ||x||_TV
%
%   .. math::  sol = arg\min_{z} \frac{1}{2} \|x - z\|_2^2 + \gamma  \|x\|_{TV}
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
%   * *param.use_fast* : Use the fast algorithm of Laurent Condat.
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
%   References: condat2013direct beck2009fastTV


% Author: Nathanael Perraudin, Marie Dankova
% Date: Jan 2013, June 2015


% Start the time counter
t1 = tic;

% for the GPU
global GLOBAL_useGPU; 

% Optional input arguments

if nargin<3, param=struct; end

if ~isfield(param, 'tol'), param.tol = 10e-4; end
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'use_fast'), param.use_fast = 0; end
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


if param.use_fast

    sol = zeros(size(x));
    
    N = size(x,1);
    for ii = 1:size(x,2)

        % a)
        k = 1; k0 = 1; kplus = 1; kminus = 1; 
        vmin = x(1,ii) - gamma; vmax = x(1,ii) + gamma;
        umin = gamma; umax = - gamma;

        while 1
            % b)
            if k == N
                sol(N,ii) = vmin + umin;
                break;
            end

            % b1)
            if x(k+1,ii) + umin < vmin - gamma
                for jj = k0:kminus
                    sol(jj,ii) = vmin;
                end
                k = kminus + 1; k0 = kminus + 1; kplus = kminus + 1; kminus = kminus + 1;
                vmin = x(k,ii);
                vmax = x(k,ii) + 2*gamma;
                umin = gamma;
                umax = - gamma;
            else
                % b2)
                if x(k+1,ii) + umax > vmax + gamma
                    for jj = k0:kplus
                        sol(jj,ii) = vmax;
                    end
                    k = kplus + 1; k0 = kplus + 1; kminus = kplus + 1; kplus = kplus + 1;
                    vmin = x(k,ii) - 2*gamma;
                    vmax = x(k,ii);
                    umin = gamma;
                    umax = - gamma;
                else
                    % b3)
                    k = k + 1;
                    umin = umin + x(k,ii) - vmin;
                    umax = umax + x(k,ii) - vmax;
                    % b31)
                    if umin >= gamma
                        vmin = vmin + (umin - gamma)/(k - k0 + 1);
                        umin = gamma;
                        kminus = k;
                    end
                    % b32)
                    if umax <= -gamma
                        vmax = vmax + (umax + gamma)/(k - k0 + 1);
                        umax = -gamma;
                        kplus = k;
                    end
                end
            end
            % c)
            if k >= N
                % c1)
                if umin < 0
                    for jj = k0:kminus
                        sol(jj,ii) = vmin;
                    end
                    k = kminus + 1; k0 = kminus +1; kminus = kminus + 1;
                    vmin = x(k,ii);
                    umin = gamma;
                    umax = x(k,ii) + gamma - vmax;
                else
                    % c2)
                    if umax > 0
                        for jj = k0:kplus
                            sol(jj,ii) = vmax;
                        end
                        k = kplus + 1;   k0 = kplus + 1; kplus = kplus + 1;
                        vmax = x(k,ii);
                        umax = - gamma;
                        umin = x(k,ii) - gamma - vmin;
                    else
                        % c3)
                        for jj = k0:N
                            sol(jj,ii) = vmin + umin/(k - k0 + 1);
                        end
                        break;
                    end
                end
            end
        end
    end

    obj = .5*norm(x(:)-sol(:), 2)^2 + gamma * sum(norm_tv1d(sol));
    crit = 'ONE-SHOT';
    rel_obj = 0;
    iter = 1;
    info.iter=iter;
    info.final_eval=obj;
    
else

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
        r = gradient_op1d(x*0);
        pold = r; 
        told = gpuArray(1); prev_obj = gpuArray(0); 
        verbose=gpuArray(param.verbose);
        tol=gpuArray(param.tol);
    else
        % Initializations
        r = gradient_op1d(x*0);
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
        sol = x - gamma*div_op1d(r);

        % Objective function value
        obj = .5*norm(x(:)-sol(:), 2)^2 + gamma * sum(norm_tv1d(sol));
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
end


if param.verbose >= 1
    if param.useGPU
        fprintf(['  GPU Prox_TV1D: obj = %e, rel_obj = %e,' ...
            ' %s, iter = %i\n'], obj, rel_obj, crit, iter);
    else
        fprintf(['  Prox_TV1D: obj = %e, rel_obj = %e,' ...
            ' %s, iter = %i\n'], obj, rel_obj, crit, iter);
    end
end




info.algo=mfilename;
info.crit=crit;
info.time=toc(t1);

end
