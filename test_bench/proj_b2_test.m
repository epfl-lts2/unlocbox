function sol = proj_b2_test(x, y, A, At, param)
% PROJ_B2 - Projection onto a L2-ball
%
% SOL = FISTA_TV_B2(xc, y, A, At, epsilon, param) solves:
% 
% 
%   min_{x} ||x - xc||_2^2   s.t.  ||y - Ax||_2 < epsilon
%
% y: measurements.
%
% A: Forward operator.
%
% At: Adjoint operator.
%
% param a Matlab structure containing the following fields:
%
%   - epsilon: Radius of the L2 ball (default = 0).
%
%   - tight: 1 if A is a tight frame or 0 if not (default = 1)
%
%   - nu: bound on the norm of the operator A, i.e.
%       ||A x||^2 <= nu * ||x||^2 (default: 1)
%
%   - tol: tolerance for the projection onto the L2 ball. The algorithms
%   stops if
%       epsilon/(1-tol) <= ||y - A z||_2 <= epsilon/(1+tol)
%   (default: 1e-3)
%
%   - max_iter: max. nb. of iterations (default: 200).
%
%   - verbose: 0 no log, 1 a summary at convergence, 2 print main
%   steps (default: 1)
%
% Author: Gilles Puy
% E-mail: gilles.puy@epfl.ch
% Date: Dec. 1, 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Optional input arguments
if ~isfield(param, 'epsilon'), param.epsilon = 0; end
if ~isfield(param, 'tight'), param.tight = 1; end
if ~isfield(param, 'tol'), param.tol = 1e-3; end
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'nu'), param.nu = 1; end
if ~isfield(param, 'max_iter'), param.max_iter = 200; end

% Useful functions for the projection
sc = @(z) z*min(param.epsilon/norm(z(:)), 1); % scaling

% Projection
if param.tight % TIGHT FRAME CASE
    
    temp = A(x) - y;
    sol = x + 1/param.nu * At(sc(temp)-temp);
    crit_B2 = 'TOL_EPS'; iter = 0;
    
else % NON TIGHT FRAME CASE
    
    % Initializations
    sol = x; u = zeros(size(y));
    iter = 1; true = 1;
    
    % Tolerance onto the L2 ball
    epsilon_low = param.epsilon/(1+param.tol);
    epsilon_up = param.epsilon/(1-param.tol);
    
    % Check if we are in the L2 ball
    norm_res = norm(y(:)-A(sol), 2);
    if norm_res <= epsilon_up
        crit_B2 = 'TOL_EPS'; true = 0;
    end
    
    % Projection onto the L2-ball
    % Init
    if param.verbose > 1
        fprintf('  Proj. B2:\n');
    end
    while true
        
        % Residual
        res = A(sol) - y; norm_res = norm(res(:), 2);
        
        % Scaling for the projection
        res = u*param.nu + res; norm_proj = norm(res(:), 2);
        
        % Log
        if param.verbose>1
            fprintf('   Iter %i, epsilon = %e, ||y - Ax||_2 = %e\n', ...
                iter, param.epsilon, norm_res);
        end
        
        % Stopping criterion
        if (norm_res>=epsilon_low && norm_res<=epsilon_up)
            crit_B2 = 'TOL_EPS'; break;
        elseif iter >= param.max_iter
            crit_B2 = 'MAX_IT'; break;
        end
        
        % Projection onto the L2 ball
        ratio = min(1, param.epsilon/norm_proj);
        u = 1/param.nu * (res - res*ratio);
        
        % Current estimate
        sol = x - At(u);
        
        % Update number of iteration
        iter = iter + 1;
        
    end
end

% Log after the projection onto the L2-ball
if param.verbose >= 1
    temp = A(sol);
    fprintf(['  Proj. B2: epsilon = %e, ||y-Ax||_2 = %e,', ...
        ' %s, iter = %i\n'], param.epsilon, norm(y(:)-temp(:)), ...
        crit_B2, iter);
end


end
