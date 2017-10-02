function [sol,info] = prox_l1(x, gamma, param)
%PROX_L1 Proximal operator with L1 norm
%   Usage:  sol=prox_l1(x, gamma)
%           sol=prox_l1(x, gamma, param)
%           [sol, info]=prox_l1(x, gamma, param)
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
%   .. sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * ||A z - y ||_1
%
%   .. math::  sol = \min_{z} \frac{1}{2} \|x - z\|_2^2 + \gamma  \|A z - y\|_1
%
%   param is a Matlab structure containing the following fields:
%
%   * *param.A* : Forward operator (default: Id).
%
%   * *param.At* : Adjoint operator (default: Id).
%
%   * *param.y* : y
%
%   * *param.tight* : 1 if A is a tight frame or 0 if not (default = 0)
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
%   We implemented the algo of  "M.J. Fadili and J-L. Starck, "Monotone
%   operator splitting for optimization problems in sparse recovery" see
%   references. See lemma 2 (section 3). The parameter nu is changed to
%   $nu^{-1}$.
%
%   See also:  proj_b1 prox_linf1 prox_l21 prox_tv
%
%   References: fadili2009monotone boyd2011distributed van2008probing

% Author: Gilles Puy, Nathanael Perraudin
% Date: Nov 2012
% Testing: test_prox_l1

% Start the time counter
t1 = tic;

% Optional input arguments
if nargin<3, param=struct; end

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'tight'), param.tight = 0; end
if ~isfield(param, 'nu'), param.nu = 1; end
if ~isfield(param, 'tol'), param.tol = 1e-3; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'At'), param.At = @(x) x; end
if ~isfield(param, 'A'), param.A = @(x) x; end
if ~isfield(param, 'weights'), param.weights = 1; end
if ~isfield(param, 'pos'), param.pos = 0; end
%    if ~isfield(param, 'y'), param.y = zeros(size(param.A(x))); end
% param.y is intitialized below to perform less computations

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




% Projection
if param.tight && ~param.pos % TIGHT FRAME CASE
    temp = param.A(x);
    if ~isfield(param, 'y'), param.y = zeros(size(temp)); end
    % TODO: improve this code to suppress the unuseful operation to compute
    % the norm_l1...
    % Done!!!
    temp2 = param.y + soft_threshold(temp -  param.y , ...
                        gamma*param.nu*param.weights) - temp;
    sol = x + 1/param.nu * param.At(temp2);
    crit = 'REL_OBJ'; iter = 1;
    dummy = temp2 +temp;
    norm_l1 = gamma*sum(param.weights(:).*abs(dummy(:)));

else % NON TIGHT FRAME CASE OR CONSTRAINT INVOLVED
    
    % Initializations
    u_l1 = zeros(size(param.A(x)));
    if ~isfield(param, 'y'), param.y = zeros(size(u_l1)); end

    sol = x - param.At(u_l1);
    prev_obj = 0; iter = 0;
    
    if param.verbose > 1
        fprintf('  Proximal l1 operator:\n');
    end
    while 1
        
        % L1 norm of the estimate
        dummy = param.A(sol);
        norm_l1 = gamma*sum(param.weights(:).*abs(dummy(:)-param.y(:)));
        norm_obj = .5*norm(x(:) - sol(:), 2)^2 + norm_l1;
        rel_obj = abs(norm_obj-prev_obj)/norm_obj;
        
        
        % Log
        if param.verbose>1
            fprintf('   Iter %i, ||A x-y||_1 = %e, rel_obj = %e\n', ...
                iter, norm_l1, rel_obj);
        end
        
        % Stopping criterion
        if (rel_obj < param.tol)
            crit = 'REL_OB'; break;
        elseif iter >= param.maxit
            crit = 'MAX_IT'; break;
        end
        
        if isnan(rel_obj)
        fprintf(['WARNING: your signal is zero. Check your timestep',...
            'or your starting point.\n']);
            crit= 'ZERO_VEC';
            break;
        end

        
        % Soft-thresholding
        res = u_l1*param.nu + dummy;
        dummy = soft_threshold(res-param.y, gamma*param.nu*param.weights)...
                + param.y;
        if param.pos
            dummy = real(dummy); dummy(dummy<0) = 0;
        end
        u_l1 = 1/param.nu * (res - dummy);
        sol = x - param.At(u_l1);
        
        
        % for comprehension of Nathanael only
        % sol=x - param.At(ul1+A(sol)-soft_threshold(ul1+A(sol)))
        
        % Update
        prev_obj = norm_obj;
        iter = iter + 1;
        
    end
    
end

% Log after the prox l1
if param.verbose >= 1
    fprintf(['  prox_L1: ||A x-y||_1 = %e,', ...
        ' %s, iter = %i\n'], norm_l1, crit, iter);
end


info.algo = mfilename;
info.iter = iter;
info.final_eval = norm_l1;
info.crit = crit;
info.time = toc(t1);
end


