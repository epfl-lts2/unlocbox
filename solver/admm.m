function [sol, info] = admm(x_0,f1, f2, param)
%ADMM alternating-direction method of multipliers
%   Usage: sol = admm(x_0,f1,f2,param);
%          sol = admm(x_0,f1,f2);
%          [sol,info,objective] = admm(...);
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
%   `admm` (using alternating-direction method of multipliers) solves:
%
%   .. sol = argmin f1(x) + f2(y) such that y=Lx
%
%   .. math::  sol = \min_x f_1(y) + f_2(x) \hspace{1cm} s.t. \hspace{1cm}  y=Lx \\
%   
%   where  $x$ is the optimization variable.
%
%   Please read the paper of Boyd "Distributed Optimization and Statistical
%   Learning via the Alternating Direction Method of Multipliers" to be
%   able to understand this demonstration file. 
%
%   *f1* is a structure representing a convex function. Inside the structure, there
%   have to be the prox of the function that can be called by *f1.proxL* and 
%   the function itself that can be called by *f1.eval*. 
%   WARNING !!!  The proxL of *f1* is not the usual prox! But the solution
%   to this problem: 
%
%   .. prox_{f1, gamma }^L(z)=min_x  1/2 ||Lx-z||_2^2 + gamma f1(x)
%
%   .. math:: prox_{f_1, \gamma }^L(z)=\min_x  \frac{1}{2} \|Lx-z\|_2^2 + \gamma f_1(x)
%
%   *f2* is a structure representing a convex function. Inside the structure, there
%   have to be the prox of the function that can be called by *f2.prox* and 
%   the function itself that can be called by *f2.eval*.
%   The prox of *f2* is the usual prox:
%
%   .. prox_{f2, gamma }(z)=min_x  1/2 ||x-z||_2^2 + gamma f2(x)
%
%   .. math:: prox_{f_2, \gamma }(z)=\min_x  \frac{1}{2} \|x-z\|_2^2 + \gamma f_2(x)           
%
%   *param* a Matlab structure containing solver paremeters. See the
%   function |solvep| for more information. Additionally it contains those
%   aditional fields:  
%
%   * *param.L* : linear operator that link $x$ and $y$: $y=Lx$. This
%     operator can be given in a matrix form (default Identity) or as a
%     function handle.
%
%   See also: solvep sdmm ppxa generalized_forward_backward
%
%   Demos:  demo_admm
%
%   References: boyd2011distributed combettes2011proximal
 
% Author: Nathanael Perraudin
% Date: 23 May 2013
% Testing: test_solvers

param.algo = 'ADMM';
[sol, info] = solvep(x_0,{f1,f2},param);

end


% function [sol, info,objective] = admm(x_0,f1, f2, param)
% %ADMM alternating-direction method of multipliers
% %   Usage: sol = admm(x_0,f1,f2,param);
% %          sol = admm(x_0,f1,f2);
% %          [sol,info,objective] = admm(...);
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
% %   `admm` (using alternating-direction method of multipliers) solves:
% %
% %   .. sol = argmin f1(x) + f2(y) such that y=Lx
% %
% %   .. math::  sol = \min_x f_1(y) + f_2(x) \hspace{1cm} s.t. \hspace{1cm}  y=Lx \\
% %   
% %   where
% %   $x$ is the variable.
% %
% %   * *x_0* is the starting point.
% %
% %   * *f1* is a structure representing a convex function. Inside the structure, there
% %     have to be the prox of the function that can be called by *f1.prox* and 
% %     the function itself that can be called by *f1.eval*. 
% %     WARNING !!!  The prox of *f1* is not the usual prox! But the solution to this problem:
% %
% %     .. prox_{f1, gamma }^L(z)=min_x  1/2 ||Lx-z||_2^2 + gamma f1(x)
% %
% %     .. math:: prox_{f_1, \gamma }^L(z)=\min_x  \frac{1}{2} \|Lx-z\|_2^2 + \gamma f_1(x)
% %
% %   * *f2* is a structure representing a convex function. Inside the structure, there
% %     have to be the prox of the function that can be called by *f2.prox* and 
% %     the function itself that can be called by *f2.eval*.
% %     The prox of *f2* is the usual prox:
% %
% %     .. prox_{f2, gamma }(z)=min_x  1/2 ||x-z||_2^2 + gamma f2(x)
% %
% %     .. math:: prox_{f_2, \gamma }(z)=\min_x  \frac{1}{2} \|x-z\|_2^2 + \gamma f_2(x)           
% %
% %   * *param* a Matlab structure containing the following fields:
% %
% %     General parameters:
% %
% %     * *param.gamma* : is the convergence parameter. By default, it's 1.
% %       (greater than 0) 
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
% %     * *param.L* : linear operator that link $x$ and $y$: $y=Lx$. This
% %       operator can be given in a matrix  form (default Identity)
% % 
% %     * *param.verbose* : 0 no log, 1 print main steps, 2 print all steps
% %       (default 1).
% %
% %     * *param.abs_tol* : If activated, this stopping critterion is the
% %       objective function being smaller than *param.tol* (default 0).
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
% %   See also:  sdmm, ppxa, generalized_forward_backward
% %
% %   Demos:  demo_admm
% %
% %   References: combettes2011proximal
% 
%  
% % Author: Nathanael Perraudin
% % Date: 23 May 2013
% % Testing: test_solver
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
% if ~isfield(param, 'gamma'), param.gamma=1 ; end
% if ~isfield(param, 'abs_tol'), param.abs_tol=1 ; end
% if ~isfield(param, 'L'), param.L=@(x) x; end
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
% % Initialization
% 
% curr_norm = f1.eval(x_0)+f2.eval(OpL(x_0));  
% [~,~,prev_norm,~,~,~] = convergence_test(curr_norm);
% [~,~,prev_rel_dual,iter,objective,~] = convergence_test(1);
% 
% y_n = OpL(x_0);
% y_old=y_n;
% z_n = zeros(size(y_n));
% 
% 
% % Main loop
% while 1
%     
%     %
%     if param.verbose >= 2
%         fprintf('Iteration %i:\n', iter);
%     end
%     
% 
%     % Algorithm
%     x_n=f1.prox(y_n-z_n,param.gamma);
%     s_n=OpL(x_n);
%     y_n=f2.prox(s_n+z_n,param.gamma);
%     reldual = norm(y_old(:) - y_n(:)) / norm(y_n(:));
% 
%     
%     z_n=z_n+s_n-y_n ;% updates
%     sol=x_n; 
%     y_old=y_n;
%     
%     % Global stopping criterion
%     curr_norm = f1.eval(sol)+f2.eval(OpL(sol));  
%     [~,rel_norm,prev_norm,~,~,~] = convergence_test(curr_norm,prev_norm);
%     [stop,~,prev_rel_dual,iter,objective,crit] = convergence_test(reldual,...
%             prev_rel_dual,iter,objective,param);
%     [sol,param] = post_process(sol, iter, curr_norm, prev_norm, objective, param);
%     if stop
%         break;
%     end
%     if param.verbose >= 2
%         fprintf(' ||f|| = %e, rel_norm = %e\n Maximum relative distance of dual variable: %e\n', ...
%             curr_norm, rel_norm, reldual);
%     end
%     
% end
% 
% % Log
% if param.verbose>=2
%     fprintf('\n Solution found:\n');
%     fprintf(' Max rel dist of dual variables: %e\n', rel_norm );
%     
%     
%     % Stopping criterion
%     fprintf(' %i iterations\n', iter);
%     fprintf(' Stopping criterion: %s \n\n', crit);
% elseif param.verbose>=1
%     fprintf('  Solution found: ||f|| = %e, rel_norm = %e, %s\n', ...
%                     curr_norm, rel_norm,crit);
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
