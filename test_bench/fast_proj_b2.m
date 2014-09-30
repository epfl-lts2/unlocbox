function [sol, infos] = fast_proj_b2(x, ~, param)
%FAST_PROJ_B2 Projection onto a L2-ball
%   Usage:  sol=fast_proj_b2(x, ~, param)
%           [sol, infos]=fast_proj_b2(x, ~, param)
%
%   Input parameters:
%         x     : Input signal.
%         param : Structure of parameters.
%   Output parameters:
%         sol   : Solution.
%         infos : Structure summarizing informations at convergence
%
%   `fast_proj_b2(x,~,param)` solves:
%
%   .. sol = argmin_{z} ||x - z||_2^2   s.t.  ||y - A z||_2 < epsilon
%
%   .. math::  sol = \min_z ||x - z||_2^2 \hspace{1cm} s.t. \hspace{1cm}  \|y - A z\|_2 < \epsilon
%
%   Remark: the projection is the proximal operator of the indicative function of
%   $\|y - A z\|_2 < \epsilon$. So it can be written:
%
%   .. prox_{f, gamma }(x)      where       f= i_c(||y - A z||_2 < epsilon)
%
%   .. math:: prox_{f, \gamma }(x) \hspace{1cm} where \hspace{1cm} f= i_c(\|y - A z\|_2 < \epsilon)
%
%
%   param is a Matlab structure containing the following fields:
%
%   * *param.y* : measurements (default: 0).
%
%   * *param.A* : Forward operator (default: Id).
%
%   * *param.At* : Adjoint operator (default: Id).
%
%   * *param.epsilon* : Radius of the L2 ball (default = 1e-3).
%
%   * *param.tight* : 1 if A is a tight frame or 0 if not (default = 1)
%
%   * *param.nu* : bound on the norm of the operator A (default: 1), i.e.
%
%   .. ` ||A x||^2 <= nu * ||x||^2 
%
%   .. math::  \|A x\|^2 \leq \nu * \|x\|^2 
%
%   * *param.tol* : tolerance for the projection onto the L2 ball  (default: 1e-3) . The algorithms
%     stops if
%   
%   .. epsilon/(1-tol) <= ||y - A z||_2 <= epsilon/(1+tol)
%
%   * *param.maxit* : max. nb. of iterations (default: 200).
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%
%   infos is a Matlab structure containing the following fields:
%
%   * *infos.algo* : Algorithm used
%
%   * *param.iter* : Number of iteration
%
%   * *param.time* : Time of execution of the function in sec.
%
%   * *param.final_eval* : Final evaluation of the function
%
%   * *param.crit* : Stopping critterion used 
%
%   * *param.residue* : Final residue  
%
%   Rem: The input "~" is useless but needed for compatibility issue.
%
%   See also:  proj_b2 proj_b1
%
%   References: fadili2009monotone beck2009fast


% Author: Gilles Puy, Perraudin Nathanael
% Date: Feb 20, 2013
%

% Start the time counter
t1 = tic;

% Optional input arguments
if ~isfield(param, 'y'), param.y = 0; end
if ~isfield(param, 'A'), param.A = @(x) x; param.At = @(x) x; end
if ~isfield(param, 'epsilon'), param.epsilon = 1e-3; end
if ~isfield(param, 'tight'), param.tight = 1; end
if ~isfield(param, 'tol'), param.tol = 1e-3; end
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'nu'), param.nu = 1; end
if ~isfield(param, 'maxit'), param.maxit = 200; end


% Useful functions for the projection
sc = @(z) z*min(param.epsilon/norm(z(:)), 1); % scaling


% Projection
if (param.tight) % TIGHT FRAME CASE   
    
    temp = param.A(x) - param.y;
    sol = x + 1/param.nu * param.At(sc(temp)-temp);
    crit = 'TOL_EPS'; iter = 0;
    u = 0;
    
else % NON TIGHT FRAME CASE
    
    % Initializations
    sol = x; u = zeros(size(param.y)); v = u;
    iter = 1; true = 1; told = 1;
    
    % Tolerance onto the L2 ball
    epsilon_low = param.epsilon/(1+param.tol);
    epsilon_up = param.epsilon/(1-param.tol);
    
    % Check if we are in the L2 ball
    dummy = param.A(sol);
    norm_res = norm(param.y(:)-dummy(:), 2);
    if norm_res <= epsilon_up
        crit = 'TOL_EPS'; true = 0;
    end
    
    % Projection onto the L2-ball
    % Init
    if param.verbose > 1
        fprintf('  Proj. B2:\n');
    end
    while true
        
        % Residual
        res = param.A(sol) - param.y; norm_res = norm(res(:), 2);
        
        % Scaling for the projection
        res = u*param.nu + res; norm_proj = norm(res(:), 2);
        
        % Log
        if param.verbose>1
            fprintf('   Iter %i, epsilon = %e, ||y - Ax||_2 = %e\n', ...
                iter, param.epsilon, norm_res);
        end
        
        % Stopping criterion
        if (norm_res>=epsilon_low && norm_res<=epsilon_up)
            crit = 'TOL_EPS'; break;
        elseif iter >= param.maxit
            crit = 'MAX_IT'; break;
        end
        
        % Projection onto the L2 ball
        t = (1+sqrt(1+4*told^2))/2;
        ratio = min(1, param.epsilon/norm_proj);
        u = v;
        v = 1/param.nu * (res - res*ratio);
        u = v + (told-1)/t * (v - u);
        
        % Current estimate
        sol = x - param.At(u);
        
        % Update number of iteration
        told = t;
        
        % Update number of iterations        
        iter = iter + 1;
        
    end
end

% Log after the projection onto the L2-ball
if param.verbose >= 1
    temp = param.A(sol);
    fprintf(['  Proj. B2: epsilon = %e, ||y-Ax||_2 = %e,', ...
        ' %s, iter = %i\n'], param.epsilon, norm(param.y(:)-temp(:)), ...
        crit, iter);
end

% Infos about algorithm
infos.algo=mfilename;
infos.iter=iter;
infos.final_eval=norm(param.A(sol) - param.y);
infos.crit=crit;
infos.residue=u;
infos.time=toc(t1);

end
