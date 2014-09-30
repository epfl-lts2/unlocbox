function [sol,info] = prox_lp(x, gamma, param)
%PROX_L1 Proximal operator with L1 norm
%   Usage:  sol=prox_lp(x, gamma)
%           sol=prox_lp(x, gamma, param)
%           [sol, info]=prox_lp(x, gamma, param)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         info  : Structure summarizing informations at convergence
%
%   `prox_l1(x, gamma, param)` solves:
%
%   .. sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * ||A z -y||_p^p
%
%   .. math::  sol = \min_{z} \frac{1}{2} \|x - z\|_2^2 + \gamma  \|A z -y \|_1^p
%
%   param is a Matlab structure containing the following fields:
%
%   * *param.p* : p norm (default: 1.5).
%
%   * *param.A* : Forward operator (default: Id).
%
%   * *param.At* : Adjoint operator (default: Id).
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
%   * *param.weights* : weights for a weighted Lp-norm (default = 1)
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
%   See also:  proj_b1 prox_l1inf prox_l12 prox_tv
%
%   References: fadili2009monotone beck2009fast


% Author: Nathanael Perraudin
% Date: Nov 2012
%

% Start the time counter
t1 = tic;

% Optional input arguments
if nargin<3, param=struct; end

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'nu'), param.nu = 1; end
if ~isfield(param, 'tol'), param.tol = 1e-3; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'At'), param.At = @(x) x; end
if ~isfield(param, 'A'), param.A = @(x) x; end
if ~isfield(param, 'y'), param.y = zeros(size(param.A(x))); end
if ~isfield(param, 'weights'), param.weights = 1; end
if ~isfield(param, 'p'), param.p = 1.5; end

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

% Useful functions
grad = @(z) z - x + param.p*gamma*param.At(sign(param.A(z)-param.y)... 
            .* abs(param.A(z)-param.y).^(param.p-1));




    
% Initializations
sol = zeros(size(x));
prev_lp = 0; iter = 0;
u_n=x;
tn=1;

stepsize=1/(param.p*gamma*max(param.weights).^2*param.nu+1);


if param.verbose > 1
    fprintf('  Proximal lp operator:\n');
end
while 1

    % Lp norm of the estimate
    dummy = param.A(sol);
    norm_lp = .5*norm(x(:) - sol(:), 2)^2 + gamma * ...
        sum(param.weights(:).*abs(dummy(:)).^param.p);
    rel_lp = abs(norm_lp-prev_lp)/norm_lp;

    % Log
    if param.verbose>1
        fprintf('   Iter %i, ||A x||_p^p = %e, rel_lp = %e\n', ...
            iter, norm_lp, rel_lp);
    end

    % Stopping criterion
    if (rel_lp < param.tol)
        crit = 'REL_OB'; break;
    elseif iter >= param.maxit
        crit = 'MAX_IT'; break;
    end

    % follow the gradient:
    
    % FISTA algorithm

    x_n=u_n-stepsize*grad(u_n);
    tn1=(1+sqrt(1+4*tn^2))/2;
    u_n=x_n+(tn-1)/tn1*(x_n-sol);
    %updates
    sol=x_n;
    tn=tn1;


    % Update
    prev_lp = norm_lp;
    iter = iter + 1;

end


% Log after the projection onto the L2-ball
if param.verbose >= 1
    fprintf(['  prox_Lp: ||A x-y||_p^p = %e,', ...
        ' %s, iter = %i\n'], norm_lp, crit, iter);
end


info.algo=mfilename;
info.iter=iter;
info.final_eval=norm_lp;
info.crit=crit;
info.time=toc(t1);
end

