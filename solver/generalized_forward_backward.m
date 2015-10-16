function [sol, info] = generalized_forward_backward(x_0, F, f , param)
%GENERALIZED_FORWARD_BACKWARD Generalized forward backward algorithm
%   Usage: sol = generalized_forward_backward(x_0,F, f2, param);
%          sol = generalized_forward_backward(x_0,F, f2);
%          [sol, info] = generalized_forward_backward(...);
%
%   Input parameters:
%         x_0   : Starting point of the algorithm
%         F     : Array of structure representing the functions to minimize
%         f2    : Another function to minimize with a known gradient
%         param : Optional parameter
%   Output parameters:
%         sol   : Solution
%         info  : Structure summarizing informations at convergence
% 
%   `generalized_forward_backward` solves:
%
%   .. sol = argmin_{z} f2(z) + Sum_{i} wi Fi(z)     for z belong to R^N
%
%   .. math::  sol = \min_{z} f_2(z) + \sum_{i} w_i F_i(z)    \hspace{1cm} for \hspace{1cm}  z\in R^N
%
%   With *z* the variable and *wi* the weight accorded to every term of the sum
%
%   * *x_0* : is the starting point.
%
%   * *F* is a cellarray of structures representing functions containing 
%     operators inside and eventually the norm. The prox: *F{i}.prox* and 
%     the function: *F{i}.eval* are defined in the same way as in the 
%     Forward-backward and Douglas-Rachford algorithms
%
%   * *f2* is a structure representing a convex function, with a  beta 
%     Lipschitz  continuous gradient. Inside the structure, there
%     have to be the gradient of the function that can be called by *f2.grad*
%     and the function itself that can be called by *f2.eval*.
%
%   * *param* is a Matlab structure containing the following fields:
%
%     * *param.weights* : weights of different functions (default = $1/N$,
%        where $N$ is the total number of function) 
%
%     * *param.lambda*: is the weight of the update term. By default 1.
%       This should be between 0 and 1.
%
%        
%   See also: solvep forward_backward ppxa sdmm
%
%   References: raguet2011generalized

% Author:  Nathanael Perraudin
% Date: Oct 24 2012
% Testing: test_solver

param.algo = 'GENERALIZED_FORWARD_BACKWARD';
if ~iscell(f)
    f = {f};
end
[sol, info] = solvep(x_0,{F{1:end},f{1:end}},param);

end


% function [sol, info,objective] = generalized_forward_backward(x_0, F, f , param)
% %GENERALIZED_FORWARD_BACKWARD Generalized forward backward algorithm
% %   Usage: sol = generalized_forward_backward(x_0,F, f2, param);
% %          sol = generalized_forward_backward(x_0,F, f2);
% %          [sol,info,objective] = generalized_forward_backward(...);
% %
% %   Input parameters:
% %         x_0   : Starting point of the algorithm
% %         F     : Array of structure representing the functions to minimize
% %         f2    : Another function to minimize with a known gradient
% %         param : Optional parameter
% %   Output parameters:
% %         sol   : Solution
% %         info  : Structure summarizing informations at convergence
% %         objective: vector (evaluation of the objectiv function each iteration)
% % 
% %   `generalized_forward_backward` solves:
% %
% %   .. sol = argmin_{z} f2(z) + Sum_{i} wi Fi(z)     for z belong to R^N
% %
% %   .. math::  sol = \min_{z} f_2(z) + \sum_{i} w_i F_i(z)    \hspace{1cm} for \hspace{1cm}  z\in R^N
% %
% %   With *z* the variable and *wi* the weight accorded to every term of the sum
% %
% %   * *x_0* : is the starting point.
% %
% %   * *F* is a cellarray of structures representing functions containing 
% %     operators inside and eventually the norm. The prox: *F{i}.prox* and 
% %     the function: *F{i}.eval* are defined in the same way as in the 
% %     Forward-backward and Douglas-Rachford algorithms
% %
% %   * *f2* is a structure representing a convex function, with a  beta 
% %     Lipschitz  continuous gradient. Inside the structure, there
% %     have to be the gradient of the function that can be called by *f2.grad*
% %     and the function itself that can be called by *f2.eval*.
% %
% %   * *param* is a Matlab structure containing the following fields:
% %
% %     General parameters:
% %
% %     * *param.gamma* : is the step size. Watch out, this parameter is bounded. It should
% %       be below $1/\beta$ (*f2* is $\beta$ Lipschitz continuous). By default, it's $10e-1$
% %
% %     * *param.weights* : weights of different functions (default = $1/N$,
% %        where $N$ is the total number of function) 
% %
% %     * *param.tol* : is stop criterion for the loop. The algorithm stops if
% %
% %       ..  (  n(t) - n(t-1) )  / n(t) < tol,
% %      
% %       .. math:: \frac{  n(t) - n(t-1) }{ n(t)} < tol,
% %
% %       where  $n(t) = f_1(Lx)+f_2(x)$ is the objective function at iteration *t*
% %       by default, `tol=10e-4`.
% %
% %     * *param.abs_tol* : If activated, this stopping critterion is the
% %       objectiv function smaller than *param.tol*. By default 0.
% %
% %     * *param.lambda*: is the weight of the update term. By default 1.
% %       This should be between 0 and 1.
% %
% %     * *param.maxit* : is the maximum number of iteration. By default, it is 200.
% % 
% %     * *param.verbose* : 0 no log, 1 print main steps, 2 print all steps.     
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
% %   See also:  douglas_rachford ppxa admm
% %
% %   Demos: demo_generalized_forward_backward
% %
% %   References: raguet2011generalized
% 
% % Author:  Nathanael Perraudin
% % Date: Oct 24 2012
% %
% 
% % Start the time counter
% t1 = tic;
% 
% % Optional input arguments
% if nargin<4, param=struct; end
% 
% if ~isfield(param, 'weights'), param.weights=ones(size(F,2),1); end
% if ~isfield(param, 'verbose'), param.verbose=1 ; end
% if ~isfield(param, 'maxit'), param.maxit=100 ; end
% if ~isfield(param, 'tol'), param.tol=1e-3 ; end
% if ~isfield(param, 'gamma'), param.gamma=1 ; end
% if ~isfield(param, 'lambda'), param.lambda=1 ; end
% 
% if nargin<3 
%     f.grad=@(x) 2*x;
%     f.eval=@(x) norm(x(:)-x_0(:),2)^2;  
% end
% 
% 
% % Normalizing the weights
% param.weights= param.weights/sum(param.weights);
% 
% 
% % Number of function
% l=size(F,2);
% 
% 
% 
% % test the evaluate function
% f = test_eval(f);
% F = test_eval(F);
% 
% 
% % Definition of the gradient function
% grad = @(y) f.grad(y);
% 
% % Algorithm - Initialisation
% z=zeros([l,size(x_0)]);
% 
% 
% for ii=1:l
%     z(ii,:,:,:,:)=x_0;
% end
% 
% 
% sol=x_0;
% curr_norm = f.eval(sol)+norm_sumg(sol,F);
% [~,~,prev_norm,iter,objective,~] = convergence_test(curr_norm);
% 
% % Algorithm - Loop
% 
% while 1
%     
%     %
%     if param.verbose > 1
%         fprintf('Iteration %i:\n', iter);
%     end
%     temp_grad=grad(sol);
%     for ii=1:length(F)
%        temp=2*sol-reshape(z(ii,:,:,:,:),size(x_0))-param.lambda*temp_grad;
%        z(ii,:,:,:,:) = reshape(z(ii,:,:,:,:),size(x_0))+ param.gamma*(F{ii}.prox(temp,1/param.weights(ii)*param.lambda)-sol);
%     end
%     
%     sol=zeros(size(x_0));
%     for ii=1:l
%         sol=sol+param.weights(ii)* reshape(z(ii,:,:,:,:),size(x_0));
%     end
%     
%     % Global stopping criterion
% 
%     curr_norm = f.eval(sol)+norm_sumg(sol,F);
%     [stop,rel_norm,prev_norm,iter,objective,crit] = convergence_test(curr_norm,prev_norm,iter,objective,param);
%     [sol, param] = post_process(sol, iter, curr_norm, prev_norm, objective, param);
%     if stop
%         break;
%     end
%     if param.verbose > 1
%         fprintf(' Norm of the general objectiv function: ||f|| = %e, rel_norm = %e\n', ...
%             curr_norm, rel_norm);
%     end
% 
% 
%   
%     
% end
% 
% 
% 
% 
% 
% % Calculation of the norm
% norm_G=f.eval(sol)+norm_sumg(sol,F);
% 
% % Log after the calculous of the prox
% if param.verbose > 1
%     fprintf('  Generalized forward backward: Sum_k ||G_k(x)|| = %e\n', norm_G);
% 
% 
%     % Stopping criterion
%     fprintf(' %i iterations\n', iter);
%     fprintf(' Stopping criterion: %s \n\n', crit);
% 
% elseif param.verbose==1
%     fprintf('  Generalized forward backward: Sum_k ||G_k(x)|| = %e', norm_G);
% 
% 
%     % Stopping criterion
%     fprintf(', it = %i', iter);
%     fprintf(', crit: %s \n', crit);
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
