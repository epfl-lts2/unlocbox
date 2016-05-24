function [sol, info] = proj_linear_eq( x,~, param )
%PROJ_LINEAR_EQ projection onto the space Az = y
%   Usage:  sol = proj_linear_eq(x, ~, param)
%           [sol, info] = proj_linear_eq(x, ~, param)
%
%   Input parameters:
%         x     : Input signal.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         info : Structure summarizing informations at convergence
%
%   `proj_linear_eq(x,~,param)` solves:
%
%   .. sol = argmin_{z} ||x - z||_2^2   s.t.  A z = y
%
%   .. math::  sol = \min_z ||x - z||_2^2 \hspace{1cm} s.t. \hspace{1cm} A z = y  
%
%   *param* is a Matlab structure containing the following fields:
%
%   * *param.y* : vector (default: 0).
%
%   * *param.method* : method used 'exact' or 'iterative' (default: 'exact').
%
%   * *param.A* : Matrix A (default: Id) (Or operator for the 'iterative'
%     method) 
%
%   * *param.At* : Matrix or operator At (Only for the 'iterative' method) 
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%   * *param.nu* : (only for iterative method) bound on the norm of the
%     operator A (default: 1), i.e.
%
%     .. ` ||A x||^2 <= nu * ||x||^2 
%
%     .. math::  \|A x\|^2 \leq \nu  \|x\|^2 
%
%   * *param.pinvA* : $A*(A A^*)^(-1)$ Pseudo inverse of A Define this
%     parameter to speed up computation (Only for 'exact').
%   
%   
%   info is a Matlab structure containing the following fields:
%
%   * *info.algo* : Algorithm used
%
%   * *info.iter* : Number of iteration
%
%   * *info.time* : Time of execution of the function in sec.
%
%   * *info.final_eval* : Final evaluation of the function
%
%   * *info.crit* : Stopping critterion used
%
%
%
%   Rem: The input "~" is useless but needed for compatibility issue.
%
%   See also:  proj_linear_ineq proj_b1
%

%
% Author: Nathanael Perraudin
% Date: May 25, 2015
% Testing: test_lp

% Start the time counter
t1 = tic;

if ~isfield(param, 'method'), param.method = 'exact'; end
if ~isfield(param, 'verbose'), param.verbose = 1; end

switch lower(param.method)
    case 'exact'
    % Optional input arguments
        if ~isfield(param, 'A'), param.A = eye(length(x)); end
        if ~isfield(param, 'pinvA'), param.pinvA = pinv(param.A); 
            warning('You should pinv A before...')
        end
        if ~isfield(param, 'verbose'), param.verbose = 1; end

        tmp = param.A*x;
        if ~isfield(param, 'y'), param.y = zeros(size(tmp)); end

        %Projection  

        sol = x - param.pinvA*(tmp-param.y);

        crit = 'TOL_EPS'; iter = 0; 
    case 'primal_dual'
        if ~isfield(param, 'A'), param.A = @(x) x; end
        if ~isfield(param, 'At'), param.At = param.A; end
        if ~isfield(param, 'nu'), param.nu = 1; end




        if isnumeric(param.A)
            A = @(x) param.A*x;
        else
            A = param.A;
        end

        if isnumeric(param.At)
            At = @(x) param.At*x;
        else
            At = param.At;
        end
        if ~isfield(param, 'y'), param.y = zeros(size(A(x))); end

        % Working fb_based_primal_dual
        f1.prox = @(z,T) param.y;
        f1.eval = @(z) eps;
        f1.L = A;
        f1.Lt = At;
        f1.norm_L = param.nu;

        f2.eval = @(z) 0.5*norm(x-z,'fro')^2;
        f2.grad = @(z) z-x;
        f2.beta = 1;

        f3.prox = @(x,T) x;
        f3.eval = @(x) 0;
        param.method = 'FISTA';

        [sol,info] = fb_based_primal_dual(x, f1,f2,f3,param);
        iter = info.iter;
        crit = info.crit;
    case 'proj_b2'
        if ~isfield(param, 'A'), param.A = @(x) x; end
        if ~isfield(param, 'At'), param.At = param.A; end




        if isnumeric(param.A)
            A = @(x) param.A*x;
        else
            A = param.A;
        end

        if isnumeric(param.At)
            At = @(x) param.At*x;
        else
            At = param.At;
        end
        if ~isfield(param, 'y'), param.y = zeros(size(A(x))); end
        paramb2 = param;
        paramb2.epsilon = 0;
        paramb2.A = A;
        paramb2.At = At;
        paramb2.tight = 0;
        paramb2.method = 'FISTA';
        [sol,info] = proj_b2(x,0,paramb2);
        iter = info.iter;
        crit = info.crit;
    otherwise
        error('Unknown method')
end

% Log after the projection
if isnumeric(param.A)
    err=norm(param.y-param.A *sol );
else
    err=norm(param.y-param.A(sol) );
end
if param.verbose >= 1
    fprintf(['  Proj. lin eq: ||y-Ax||_2 = %e,', ...
        ' %s, iter = %i\n'],err , crit, iter);
end

% info about algorithm
info.algo=mfilename;
info.iter=iter;
info.final_eval=err;
info.crit=crit;
info.time=toc(t1);

end


