function [sol, info] = chambolle_pock(x_0,f1, f2, param)
%CHAMBOLLE_POCK A First-Order Primal-Dual Algorithm by Chambolle and Pock
%   Usage: sol = chambolle_pock(x_0,f1,f2,param);
%          sol = chambolle_pock(x_0,f1,f2);
%          [sol, info] = chambolle_pock(...);
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
%   `chambolle_pock` solves:
%
%   ..  min max  <Lx;y> f1(x) - f2*(y)
%   ..   x   y 
%
%   .. math::  \min_x \max_y \langle Lx ; y \rangle + f_1(x) - f_2^*(y)  
%
%   where $x,y$ are the optimization variables and $f_2^*$ is the the
%   convex conjugate of $f_2$.
%
%   This is the dual problem of:
%
%   ..  min  f2(Lx) +  f1(x)
%   ..   x  
%
%   .. math:: \min_x f_2(Lx) + f_1(x)  
%
%   Note that *f2* is a structure of a functions with:
%
%   * `f2.eval(x_i)` : an operator to evaluate the function
%   * `f2.prox(x_i, gamma)` : an operator to evaluate the prox of the function
%
%   Optionally you can define
%
%   * `f2.L`  : linear operator, matrix or operator (default identity)
%   * `f2.Lt` : adjoint of linear operator, matrix or operator (default identity)
%   * `f2.norm_L` : bound on the norm of the operator L (default: 1), i.e.
%
%     .. ` ||L x||^2 <= nu * ||x||^2 
%
%     .. math::  \|L x\|^2 \leq \nu \|x\|^2
%
%   *param.tau* and *param.sigma* need to satisfy:
%
%   .. tau * beta * nu < 1
%
%   .. math:: \tau * \beta * \nu < 1
%
%   The default choice for the time-step makes the following 
%
%     .. tau  = beta = 1/sqrt(nu)
%
%     .. math::  {\tau} - \sigma  = \frac{1}{\sqrt{\nu}}   
%
%   *param* a Matlab structure containing solver paremeters. See the
%   function |solvep| for more information.
%
%   See also: fb_based_primal_dual, admm, sdmm
%
%   References: chambolle2010first

 
% Author: Nathanael Perraudin
% Date: 23 May 2013
% Testing: test_solvers

param.algo = 'CHAMBOLLE_POCK';
[sol, info] = solvep(x_0,{f1,f2},param);



end
 

% function [sol, info,objective] = chambolle_pock(x_0,f1, f2, param)
% %CHAMBOLLE_POCK A First-Order Primal-Dual Algorithm by Chambolle and Pock
% %   Usage: sol = chambolle_pock(x_0,f1,f2,param);
% %          sol = chambolle_pock(x_0,f1,f2);
% %          [sol,info,objective] = chambolle_pock(...);
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
% %   `chambolle_pock` solves:
% %
% %   ..  min max  <Lx;y> f1(x) - f2*(y)
% %   ..   x   y 
% %
% %   .. math::  \min_x \max_y \langle Lx ; y \rangle + f_1(x) - f_2^*(y)  
% %
% %   where $x,y$ are the optimization variables and $f_2^*$ is the the
% %   convex conjugate of $f_2$.
% %
% %   This is the dual problem of:
% %
% %   ..  min  f2(Lx) +  f1(x)
% %   ..   x  
% %
% %   .. math:: \min_x f_2(Lx) + f_1(x)  
% %
% %   The algorithm only returns the optimal variable $x$ in the field `sol`
% %
% %   * *x_0* is the starting point.
% %
% %   * *f1* is a structure representing a convex function. Inside the structure, there
% %     have to be the prox of the function that can be called by *f1.prox* and 
% %     the function itself that can be called by *f1.eval*.
% %
% %   * *f2* is a structure representing a convex function. Inside the structure, there
% %     have to be the prox of the function that can be called by *f2.prox* and 
% %     the function itself that can be called by *f2.eval*.
% %
% %   * *param* a Matlab structure containing the following fields:
% %
% %     * *param.gamma* : is the stepsize: [0,1]. By default, it's 1.
% %
% %     * *param.tau* : convergence parameter of $f_1$. It is a strictly
% %       positiv number. Be default, it's 1.
% %
% %     * *param.rho* : convergence parameter of $f_2$. It is a strictly
% %       positiv number. Be default, it's 1.
% %
% %     * *param.tol* : is stop criterion for the loop. The algorithm stops if
% %
% %       ..  (||  y(t) - y(t-1) ||)  /  || y(t) || < tol,
% %      
% %       .. math:: \frac{ \| y(t) - y(t-1) \| }{\| n(t)\|} < tol,
% %
% %       where  $y(t)$ are the dual the objective function at iteration *t*
% %       by default, `tol=10e-4`.
% %
% %     * *param.maxit* : is the maximum number of iteration. By default, it is 200.
% %
% %     * *param.L* : linear operator. This operator can be given in a matrix
% %       form (default Identity) 
% % 
% %     * *param.Lt* : adjoint operator of *param.L* (default Identity)
% %
% %     * *param.verbose* : 0 no log, 1 print main steps, 2 print all steps
% %       (default 1).
% %
% %     * *param.abs_tol* : If activated, this stopping critterion is the
% %       objectiv function being smaller than *param.tol* (default 0).
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
% %   * *param.rel_norm* : Relative norm at convergence 
% %
% %
% %   See also: admm, sdmm
% %
% %   Demos:  demo_chambolle_pock.
% %
% %   References: chambolle2010first
% 
%  
% % Author: Nathanael Perraudin
% % Date: 23 May 2013
% %
% 
% % Start the time counter
% t1 = tic;
% 
% 
% % Optional input arguments
% if nargin<4, param=struct; end
% 
% if ~isfield(param, 'tol'), param.tol=10e-4 ; end
% if ~isfield(param, 'maxit'), param.maxit=200; end
% if ~isfield(param, 'verbose'), param.verbose=1 ; end
% if ~isfield(param, 'lambda'), param.lambda=1 ; end
% if ~isfield(param, 'gamma'), param.gamma=1 ; end
% if ~isfield(param, 'tau'), param.tau=1 ; end
% if ~isfield(param, 'rho'), param.rho=1 ; end
% if ~isfield(param, 'abs_tol'), param.abs_tol=1 ; end
% if ~isfield(param, 'L'), param.L=@(x) x; end
% if ~isfield(param, 'Lt'), param.Lt=@(x) x; end
% 
% 
% 
% % test the evaluate function
% [f1] = test_eval(f1);
% [f2] = test_eval(f2);
% 
% if isa(param.L,'numeric')
%    OpL= @(x) param.L*x;
% else
%    OpL= param.L;
% end
% 
% if isa(param.Lt,'numeric')
%    OpLt= @(x) param.Lt*x;
% else
%    OpLt= param.Lt;
% end
% 
% % Initialization
% 
% curr_norm = f1.eval(x_0)+f2.eval(OpL(x_0));  
% [~,~,prev_norm,~,~,~] = convergence_test(curr_norm);
% [~,~,prev_rel_dual,iter,objective,~] = convergence_test(1);
% 
% 
% y_n = OpL(x_0);
% z_n = x_0;
% x_n=x_0;
% 
% 
% % Main loop
% while 1
%     
%     %
%     if param.verbose >= 1
%         fprintf('Iteration %i:\n', iter);
%     end
%     
% 
%     % Algorithm
%     y_n_old=y_n;
%     y_n=prox_adjoint(y_n+param.rho*OpL(z_n),param.rho,f2);
%     x_n_old=x_n;
%     x_n=f1.prox(x_n+param.tau*OpLt(y_n),param.tau);
%     z_n=x_n+param.gamma*(x_n-x_n_old);
%     reldual=norm(y_n_old-y_n)/norm(y_n);
% 
%     
%     sol=z_n; 
% 
%     
%     % Global stopping criterion
%     curr_norm = f1.eval(sol)+f2.eval(OpL(sol));  
%     [~,rel_norm,prev_norm,~,~,~] = convergence_test(curr_norm,prev_norm);
%     [stop,~,prev_rel_dual,iter,objective,crit] = convergence_test(reldual,...
%             prev_rel_dual,iter,objective,param);
%     [z_n,param] = post_process(sol, iter, curr_norm, prev_norm, objective, param);
%     if stop
%         break;
%     end
%     if param.verbose >= 1
%         fprintf(' ||f|| = %e, rel_norm = %e\n Maximum relative distance of dual variable: %e\n', ...
%             curr_norm, rel_norm, reldual);
%     end
%     
% end
% 
% % Log
% if param.verbose>=1
%     fprintf('\n Solution found:\n');
%     fprintf(' Final relative norm: %e\n', rel_norm );
%     
%     
%     % Stopping criterion
%     fprintf(' %i iterations\n', iter);
%     fprintf(' Stopping criterion: %s \n\n', crit);
%     
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

