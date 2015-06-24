function [ sol ] = proj_linear_eq( x,~, param )
%PROJ_LINEAR_EQ projection onto the space Az = y
%   Usage:  sol = proj_linear_eq(x, ~, param)
%           [sol, infos] = proj_linear_eq(x, ~, param)
%
%   Input parameters:
%         x     : Input signal.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         infos : Structure summarizing informations at convergence
%
%   `proj_dual(x,~,param)` solves:
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
%   infos is a Matlab structure containing the following fields:
%
%   * *infos.algo* : Algorithm used
%
%   * *infos.iter* : Number of iteration
%
%   * *infos.time* : Time of execution of the function in sec.
%
%   * *infos.final_eval* : Final evaluation of the function
%
%   * *infos.crit* : Stopping critterion used
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

if strcmp(param.method,'exact')
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
else
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

    % Working fb_based_primal_dual
    f1.prox = @(z,T) param.y;
    f1.eval = @(z) eps;
    f1.L = A;
    f1.Lt = At;
    
    f2.eval = @(z) 0.5*norm(x-z,'fro')^2;
    f2.grad = @(z) z-x;
    f2.beta = 1;
    
    f3.prox = @(x,T) x;
    f3.eval = 0;
    param.method = 'FISTA';
    
    [sol,infos] = fb_based_primal_dual(x, f1,f2,f3,param);

%     paramb2 = param;
%     paramb2.epsilon = 0;
%     paramb2.A = A;
%     paramb2.At = At;
%     paramb2.method = 'FISTA';
%     [sol,infos] = proj_b2(x,0,paramb2);
    iter = infos.iter;
    crit = infos.crit;
end

% Log after the projection
error=norm(param.y-param.A *sol );
if param.verbose >= 1
    fprintf(['  Proj. lin eq: ||y-Ax||_2 = %e,', ...
        ' %s, iter = %i\n'],error , crit, iter);
end

% Infos about algorithm
infos.algo=mfilename;
infos.iter=iter;
infos.final_eval=error;
infos.crit=crit;
infos.time=toc(t1);

end


