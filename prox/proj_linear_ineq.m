function [ sol ] = proj_linear_ineq( x,~, param )
%PROJ_LINEAR_INEQ projection onto the space Az = y
%   Usage:  sol = proj_linear_ineq(x, ~, param)
%           [sol, infos] = proj_linear_ineq(x, ~, param)
%
%   Input parameters:
%         x     : Input signal.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         infos : Structure summarizing informations at convergence
%
%   `proj_linear_ineq(x,~,param)` solves:
%
%   .. sol = argmin_{z} ||x - z||_2^2   s.t.  A z <= y
%
%   .. math::  sol = \min_z ||x - z||_2^2 \hspace{1cm} s.t. \hspace{1cm} A z \leq y  
%
%   param is a Matlab structure containing the following fields:
%
%   * *param.y* : measurements (default: 0).
%
%   * *param.A* : Matrix A (default: Id).
%
%   * *param.AAtinv* : $(A A^*)^(-1)$ Define this parameter to speed up
%     computation. 
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
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
%   * *infos.residue* : Final residue  
%
%
%   Rem: The input "~" is useless but needed for compatibility issue.
%
%   See also:  prox_l2 proj_b1
%
%   References: 

%
% Author: Nathanael Perraudin
% Date: May 25, 2015
% Testing: test_lp

% Start the time counter
t1 = tic;

% Optional input arguments

if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'method'), param.method = 'quadprog'; end

if ~isfield(param, 'A'), param.A = eye(length(x)); end

if isnumeric(param.A)
    A = @(x) param.A*x;
else
    A = param.A;
end
tmp = A(x);

if ~isfield(param, 'y'), param.y = zeros(size(tmp)); end

if sum(tmp-param.y > 0)
    
    if strcmp(param.method,'quadprog')
        H = eye(length(x));
        f = -x;
        A = param.A;
        b = param.y;
        [sol] = quadprog(H,f,A,b);
    elseif strcmp(param.method,'iterative')
        warning('The help for this method has not been done')
        if ~isfield(param, 'At'), param.At = param.A; end
        if ~isfield(param, 'y'), param.y = zeros(size(tmp)); end
        
        if isnumeric(param.At)
            At = @(x) param.At*x;
        else
            At = param.At;
        end

        f1.prox = @(z,T) proj_inf(z,param.y);
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
    else
        error('Unknown method')
    end
else
    sol = x;
end

 % TODO to be improved   
crit = 'TOL_EPS'; iter = 0;


% Log after the projection onto the L2-ball
error_num=norm(param.y-param.A *sol )^2;
if param.verbose >= 1
    fprintf(['  Proj. lin ineq: ||y-Ax||_2^2 = %e,', ...
        ' %s, iter = %i\n'],error_num , crit, iter);
end

% Infos about algorithm
infos.algo=mfilename;
infos.iter=iter;
infos.final_eval=error_num;
infos.crit=crit;
infos.time=toc(t1);

end


function x = proj_inf(x,y)
    x(x>y) = y(x>y);
end

