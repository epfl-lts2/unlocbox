function [sol,info] = prox_l2(x, gamma, param)
%PROX_L2 Proximal operator with L2 norm
%   Usage:  sol=prox_l2(x, gamma)
%           sol=prox_l2(x, gamma, param)
%           [sol, info]=prox_l2(x, gamma, param)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         info  : Structure summarizing informations at convergence
%
%   `prox_l2(x, gamma, param)` solves:
%
%   .. sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * ||w (A z - y)||_2^2
%
%   .. math::  sol = \min_{z} \frac{1}{2} \|x - z\|_2^2 + \gamma  \|w(A z-y)\|_2^2
%
%   where w are some weights.
%
%   param is a Matlab structure containing the following fields:
%
%   * *param.weights* : weights for a weighted L2-norm (default = 1)
%
%   * *param.y* : measurements (default: 0).
%
%   * *param.A* : Forward operator (default: Id).
%
%   * *param.At* : Adjoint operator (default: A).
%
%   * *param.tightT* : 1 if $A^T$ is a tight frame or 0 if not (default = 0)
%     Note that $A^T$ tight means $A A^T = \nu I$.
%
%   * *param.tight* : 1 if A is a tight frame or 0 if not (default = 0)
%     Note that $A$ tight means $A^T A = \nu I$.
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
%   * *param.pcg* : Use the fast PCG algorithm (default 1).
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
%   See also:  proj_b2 prox_l1
 
% Author: Nathanael Perraudin
% Date: Nov 2012
% Testing: test_prox_l2
 
% Start the time counter
t1 = tic;
 
% Optional input arguments
if nargin<3, param=struct; end
 
% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'tight'), param.tight = 0; end
if ~isfield(param, 'tightT'), param.tightT = 0; end
if ~isfield(param, 'nu'), param.nu = 1; end
if ~isfield(param, 'tol'), param.tol = 1e-3; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'A'), param.A = @(x) x; end
if ~isfield(param, 'At'), param.At = param.A; end
if ~isfield(param, 'y'), param.y = zeros(size(param.A(x))); end
if ~isfield(param, 'weights'), param.weights = ones(size(param.y)); end
if ~isfield(param, 'pcg'), param.pcg = 1; end
 
 
 
 
% Check if the weight are correct with respect to the tight frame
if (max(param.weights(:))~=min(param.weights(:))) && param.tight
    param.tight=0;    
end
 
 
if max(param.weights(:))==min(param.weights(:))
    param.weights=max(param.weights(:));
end
 
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
 
param.weights=test_weights(param.weights);
 
 
 
% Projection
if param.tight % TIGHT FRAME CASE
    
    sol=(x+gamma*2*param.At(param.y.*abs(param.weights).^2))./(gamma*2*param.nu*param.weights.^2+1);
    
    
    % Infos for log...
    % L2 norm of the estimate
    dummy = (param.A(sol)-param.y);
    norm_l2 = .5*norm(x(:) - sol(:), 'fro')^2 + gamma *  norm(param.weights(:).*dummy(:))^2;
    % stopping criterion
    crit = 'REL_OB'; 
    % number of iteration
    iter=0;
    
elseif param.tightT % TIGHT FRAME CASE OF A^T
    
    sol=(x - 2*gamma./((2*gamma*param.weights.^2 +1)*param.nu) .*...
        param.At((param.A(x)-param.y).*abs(param.weights).^2));

    
    
    % Infos for log...
    % L2 norm of the estimate
    dummy = (param.A(sol)-param.y);
    norm_l2 = .5*norm(x(:) - sol(:), 'fro')^2 + gamma *  norm(param.weights(:).*dummy(:))^2;
    % stopping criterion
    crit = 'REL_OB'; 
    % number of iteration
    iter=0;
else % NON TIGHT FRAME
      
    if param.pcg
        
        %(2*gamma*F*diag(w.^2)*F'+eye(L))^(-1)*(x+2*gamma.*F*(w.^2.*y));
%         A = @(z) 2*gamma*param.At(param.weights.^2.*param.A(z))+z;
%         b = 2*gamma*param.At(param.weights.^2.*param.y) + x;        
%         [sol,flag,~,iter] = pcg(A,b,param.tol,param.maxit,[],[],x);       
        sx = size(x);
 
        A = @(z) reshape( 2*gamma*param.At( param.weights.^2 .* ...
            param.A(reshape(z,sx))),[],1) + z;
        b = 2*gamma*param.At(param.weights.^2.*param.y) + x;        
        [sol,flag,~,iter] = pcg(A,b(:),param.tol,param.maxit,[],[],x(:));
        
        sol = reshape(sol,sx);
       
        dummy = param.weights.*(param.A(sol)-param.y);
        norm_l2 = .5*norm(x(:) - sol(:), 'fro')^2 + gamma * norm(dummy(:))^2;
        if ~flag
            crit = 'REL_OB'; 
        else 
            crit = 'MAX_IT'; 
        end
    else
        % Initializations
        u_n=x;
        sol=x;
        tn=1;
        prev_l2 = 0; iter = 0;
        % stepsize
        stepsize=1/(2*gamma*max(abs(param.weights)).^2*param.nu+1);
        % gradient
        grad= @(z) z-x+gamma*2.*param.At(param.weights.^2.*(param.A(z)-param.y));
 
        % Init
        if param.verbose > 1
            fprintf('  Proximal l2 operator:\n');
        end
        while 1
 
            % L2 norm of the estimate
            dummy = param.weights.*(param.A(sol)-param.y);
            norm_l2 = .5*norm(x(:) - sol(:), 'fro')^2 + gamma * norm(dummy(:))^2;
            rel_l2 = abs(norm_l2-prev_l2)/norm_l2;
 
            % Log
            if param.verbose>1
                fprintf('   Iter %i, ||w (A x- y)||_2^2 = %e, rel_l2 = %e\n', ...
                    iter, norm_l2, rel_l2);
            end
 
            % Stopping criterion
            if (rel_l2 < param.tol)
                crit = 'REL_OB'; break;
            elseif iter >= param.maxit
                crit = 'MAX_IT'; break;
            end
 
            % FISTA algorithm
            x_n=u_n-stepsize*grad(u_n);
            tn1=(1+sqrt(1+4*tn^2))/2;
            u_n=x_n+(tn-1)/tn1*(x_n-sol);
            %updates
            sol=x_n;
            tn=tn1;
 
 
            % Update
            prev_l2 = norm_l2;
            iter = iter + 1;
 
        end
    end
end
 
% Log after the projection onto the L2-ball
if param.verbose >= 1
    fprintf('  prox_L2: ||w (A x- y) ||_2^2 = %e, %s, iter = %i\n', ...
    norm_l2, crit, iter);
 
end
 
 
info.algo=mfilename;
info.iter=iter;
info.final_eval=norm_l2;
info.crit=crit;
info.time=toc(t1);
 
end
