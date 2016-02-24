function [sol, info] = forward_backward(x_0,f1, f2, param)
%FORWARD_BACKWARD Forward-backward splitting algorithm
%   Usage: sol = forward_backward(x_0,f1, f2, param);
%          sol = forward_backward(x_0,f1, f2);
%          [sol,infos] = forward_backward(...);
%
%   Input parameters:
%         x_0   : Starting point of the algorithm
%         f1    : First function to minimize
%         f2    : Second function to minimize
%         param : Optional parameter
%   Output parameters:
%         sol   : Solution
%         info  : Structure summarizing informations at convergence
%
%   `forward_backward` solves:
%
%   .. sol = argmin f1(x) + f2(x)      for x belong to R^N
%
%   .. math::  sol = arg \min_x f_1(x) + f_2(x) \hspace{1cm} for \hspace{1cm}  x\in R^N
%
%   where *x* is the optimization variable.
%
%   *f1* is a structure representing a convex function. Inside the structure, there
%   have to be the prox of the function that can be called by *f1.prox* and 
%   the function itself that can be called by *f1.eval*. 
%
%   *f2* is a structure representing a convex function, with a  $\beta$ 
%   Lipschitz  continuous gradient. Inside the structure, there
%   have to be the gradient of the function that can be called by *f2.grad* 
%   and the function itself that can be called by *f2.eval*.
%
%   *param* a Matlab structure containing solver paremeters. See the
%   function |solvep| for more information. Additionally it contains those
%   aditional fields:  
%
%   * *param.lambda* : is the weight of the update term. It is kind of a
%     timestep for the proximal operators. (Warning it should not be
%     confused with *gamma*, the time step for gradient descent part). By
%     default it is set to 1. Do not change this parameter unless you know
%     what you do.
%
%   * *param.method* : is the method used to solve the problem. It can be
%     the fast version 'FISTA' or 'ISTA'. By default, it's 'FISTA'. 
%
%
%   See also:  solvep douglas_rachford admm generalized_forward_backward
%
%   References: beck2009fast combettes2007douglas
%

% Author: Nathanael Perraudin
% Date: 24 oct 2012
% Testing: test_solver

param.algo = 'FORWARD_BACKWARD';
[sol, info] = solvep(x_0,{f2,f1},param);

end

% function [sol, info,objective] = forward_backward(x_0,f1, f2, param)
% %FORWARD_BACKWARD Forward-backward splitting algorithm
% %   Usage: sol = forward_backward(x_0,f1, f2, param);
% %          sol = forward_backward(x_0,f1, f2);
% %          [sol,infos,objectiv] = forward_backward(...);
% %
% %   Input parameters:
% %         x_0   : Starting point of the algorithm
% %         f1    : First function to minimize
% %         f2    : Second function to minimize
% %         param : Optional parameter
% %   Output parameters:
% %         sol   : Solution
% %         info  : Structure summarizing informations at convergence
% %         objective: vector (evaluation of the objectiv function each iteration)
% %
% %   `forward_backward` solves:
% %
% %   .. sol = argmin f1(x) + f2(x)      for x belong to R^N
% %
% %   .. math::  sol = arg \min_x f_1(x) + f_2(x) \hspace{1cm} for \hspace{1cm}  x\in R^N
% %
% %   where *x* is the variable.
% %
% %   * *x_0* is the starting point.
% %
% %   * *f1* is a structure representing a convex function. Inside the structure, there
% %     have to be the prox of the function that can be called by *f1.prox* and 
% %     the function itself that can be called by *f1.eval*. 
% %
% %   * *f2* is a structure representing a convex function, with a  $\beta$ 
% %     Lipschitz  continuous gradient. Inside the structure, there
% %     have to be the gradient of the function that can be called by *f2.grad* 
% %     and the function itself that can be called by *f2.eval*.
% %
% %   * *param* a Matlab structure containing the following fields:
% %
% %     General parameters:
% %
% %     * *param.gamma* : is the step size. Watch out, this parameter is bounded. It should
% %       be below $1/\beta$ (*f2* is $\beta$ Lipchitz continuous). By default, it's $10e-1$
% %
% %     * *param.tol* : is stop criterion for the loop. The algorithm stops if
% %
% %       ..  (  n(t) - n(t-1) )  / n(t) < tol,
% %      
% %       .. math:: \frac{  n(t) - n(t-1) }{ n(t)} < tol,
% %
% %       where  $n(t) = f_1(x)+f_2(x)$ is the objective function at iteration *t*
% %       by default, `tol=10e-4`.
% %
% %     * *param.method* : is the method used to solve the problem. It can be 'FISTA' or
% %       'ISTA'. By default, it's 'FISTA'.
% %
% %     * *param.lambda*: is the weight of the update term in ISTA method. By default 1.
% %       This should be between 0 and 1. It's set to 1 for FISTA.
% %
% %     * *param.maxit* : is the maximum number of iteration. By default, it is 200.
% % 
% %     * *param.verbose* : 0 no log, 1 print main steps, 2 print all steps.
% %
% %     * *param.abs_tol* : If activated, this stopping critterion is the
% %       objectiv function smaller than *param.tol*. By default 0.
% %
% %
% %   info is a Matlab structure containing the following fields:
% %
% %   * *info.algo* : Algorithm used
% %
% %   * *info.iter* : Number of iteration
% %
% %   * *info.time* : Time of exectution of the function in sec.
% %
% %   * *info.final_eval* : Final evaluation of the objectivs functions
% %
% %   * *info.crit* : Stopping critterion used 
% %
% %   * *info.rel_norm* : Relative norm at convergence 
% %
% %             
% %   See also:  douglas_rachford admm generalized_forward_backward
% %
% %   Demos: demo_forward_backward 
% %
% %   References: beck2009fast combettes2007douglas
% 
% 
% % Author: Nathanael Perraudin
% % Date: 24 oct 2012
% 
% % Start the time counter
% t1 = tic;
% 
% % Optional input arguments
% if nargin<4, param=struct; end
% 
% if ~isfield(param, 'gamma'), param.gamma = 1; end
% if ~isfield(param, 'tol'), param.tol=10e-4 ; end
% if ~isfield(param, 'maxit'), param.maxit=200; end
% if ~isfield(param, 'verbose'), param.verbose=1 ; end
% if ~isfield(param, 'lambda'), param.lambda=1 ; end
% if ~isfield(param, 'method'), param.method='FISTA' ; end
% 
% if nargin<3 
%     f2.grad=@(x) 2*(x-x_0);
%     f2.eval=@(x) norm(x(:)-x_0(:),2)^2;  
% end
% 
% if nargin<2
%     f1.prox=@(x) prox_L1(x, 1, param);
%     f1.eval=@(x) norm(x(:),1);   
% end;
% 
% % test the evaluate function
% [f1] = test_eval(f1);
% [f2] = test_eval(f2);
% 
% 
% 
% % Initialization
%     curr_norm = f1.eval(x_0)+f2.eval(x_0);  
% [~,~,prev_norm,iter,objective,~] = convergence_test(curr_norm);
% 
% if param.verbose>=2, 
%     if strcmp(param.method,'ISTA'),
%         fprintf('Algorithm selected: ISTA\n');
%         
%     else
%         fprintf('Algorithm selected: FISTA\n');
% 
%     end
% end
% % ISTA
% x_n = x_0;
% 
% %FISTA
% u_n=x_0;
% sol=x_0;
% tn=1;
% 
% % Main loop
% while 1
%     
%     %
%     if param.verbose >= 2
%         fprintf('Iteration %i:\n', iter);
%     end
%     
%     if strcmp(param.method,'ISTA'),
%         % ISTA algorithm
%         u_n=x_n-param.gamma*f2.grad(x_n);
%         x_n=x_n+param.lambda*(f1.prox(u_n, param.gamma)-x_n);
%         sol = x_n; % updates
%     else
%         % FISTA algorithm
%         x_n=f1.prox(u_n-param.gamma*f2.grad(u_n), param.gamma);
%         tn1=(1+sqrt(1+4*tn^2))/2;
%         u_n=x_n+(tn-1)/tn1*(x_n-sol);
%         %updates
%         sol=x_n;
%         tn=tn1;
%     end
%     
%     % Global stopping criterion
%     curr_norm = f1.eval(sol)+f2.eval(sol);  
%     [stop,rel_norm,prev_norm,iter,objective,crit] = convergence_test(curr_norm,prev_norm,iter,objective,param);
%     [x_n,param] = post_process(sol, iter, curr_norm, prev_norm, objective, param);
%     if param.verbose >= 2
%         fprintf('  ||f|| = %e, rel_norm = %e\n', ...
%             curr_norm, rel_norm);
%     end
%     if stop
%         break;
%     end
%     
% end
% 
% % Log
% if param.verbose>=2
%     % Print norm
%     fprintf('\n Forward backward:\n');
%     fprintf(' Final relative norm: %e\n', rel_norm );    
%     % Stopping criterion
%     fprintf(' %i iterations\n', iter);
%     fprintf(' Stopping criterion: %s \n\n', crit);
% elseif param.verbose>=1
%     fprintf('  Forward backward: ||f|| = %e, rel_norm = %e, it = %i, %s\n', ...
%                     curr_norm, rel_norm, iter,crit);
% end
% 
% info.algo=mfilename;
% info.iter=iter;
% info.final_eval=curr_norm;
% info.crit=crit;
% info.time=toc(t1);
% info.rel_norm=rel_norm;
% 
% end
