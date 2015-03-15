function [sol, info] = prox_sumg(x, gamma , param)
%PROX_sumG Proximal operator of a sum of function
%   Usage:  sol=prox_sumg(x, gamma, param)
%           [sol, info]=prox_sumg(...)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of parameters.
%   Output parameters:
%         sol   : Solution.
%         info  : Structure summarizing informations at convergence
%
%   `prox_sumG(x, gamma , param)` solves:
%
%   .. sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * Sum_{i} w_i G_i(z)     for z,x belong to R^N
%
%   .. math::  sol = \operatorname*{arg\,min}_{z} \frac{1}{2} \|x - z\|_2^2 + \gamma  \sum_{i} w_i G_i(z)      \hspace{1cm} for \hspace{1cm}  z,x\in R^N
%
%   param is a Matlab structure containing the following fields:
%
%   * *param.G* : cellarray of structure with all the prox operator inside and eventually 
%     the norm if no norm is defined, the $L^1$ norm is used the prox: *F{i}.prox* and 
%     norm: *F{i}.eval* are defined in the same way as in the
%     Forward-backward and Douglas-Rachford algorithms
%
%   * *param.weights* : weights of different functions (default = $1/N$,
%     where $N$ is the total number of function) 
%
%   * *param.tol* : is stop criterion for the loop. The algorithm stops if
%
%     ..  (  n(t) - n(t-1) )  / n(t) < tol,
%      
%     .. math:: \frac{  n(t) - n(t-1) }{ n(t)} < tol,
%
%     where  $n(t) = f_1(Lx)+f_2(x)$ is the objective function at iteration *t*
%     by default, `tol=10e-4`.
%
%   * *param.lambda_t*: is the weight of the update term. By default 1.
%     This should be between 0 and 1.
%
%   * *param.maxit* : is the maximum number of iteration. By default, it is 200.
% 
%   * *param.verbose* : 0 no log, 1 print main steps, 2 print all steps.     
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
%   See also:  generalized_forward_backward
%
%   Demo: demo_prox_multi_functions
%
%   References: raguet2011generalized


% Author:  Nathanael Perraudin
% Date: October 2012
%


% Start the time counter
t1 = tic;

% Optional input arguments
if nargin<3, param=struct; end

if ~isfield(param, 'G'), param.G = struct([]); end
if ~isfield(param, 'weights'), param.weights=ones(length(param.G),1); end
if ~isfield(param, 'verbose'), param.verbose=1 ; end
if ~isfield(param, 'maxit'), param.maxit=100 ; end
if ~isfield(param, 'tol'), param.tol=1e-3 ; end

% gamma should maybe be always set to 1 is this function
param.gamma=1;
if ~isfield(param, 'lambda_t'), param.lambda_t=1 ; end

% Test of gamma
if test_gamma(gamma)
    sol = x;
    info.algo=mfilename;
    info.iter=0;
    info.final_eval=0;
    info.crit='--';
    info.time=toc(t1);
    return; 
end

% Test of the weights
param.weights=test_weights(param.weights);

% Normalizing the weights
param.weights= param.weights/sum(param.weights);


% Number of function
l=length(param.G);




% Definition of the gradient function
grad = @(y) (y-x);

% Algorithm - Initialisation
z=zeros([l,size(x)]);


for ii=1:l
    z(ii,:,:)=x;
end


sol=x;
iter=1;

prev_norm = 0.5*norm(sol,2)^2+gamma*norm_sumg(sol,param.G);

% Algorithm - Loop

while 1
    
    %
    if param.verbose >= 1
        fprintf('Iteration %i:\n', iter);
    end
    
    for ii=1:l
       temp=2*sol-reshape(z(ii,:,:),size(x))-param.lambda_t*grad(sol);
       z(ii,:,:) = reshape(z(ii,:,:),size(x))+ param.gamma* ...
       (param.G{ii}.prox(temp,gamma/param.weights(ii)*param.lambda_t)-sol);
    end
    
    sol=zeros(size(x));
    for ii=1:l
        sol=sol+param.weights(ii)* reshape(z(ii,:,:),size(x));
    end
    
    % Global stopping criterion

    curr_norm = 0.5*norm(sol,2)^2+gamma*norm_sumg(sol,param.G);
    rel_norm = abs(curr_norm - prev_norm)/curr_norm;
    if param.verbose >= 1
        fprintf('  ||f|| = %e, rel_norm = %e\n', ...
            curr_norm, rel_norm);
    end
    if (rel_norm < param.tol)
        crit = 'REL_NORM';
        break;
    elseif iter >= param.maxit
        crit = 'MAX_IT';
        break;
    end
    
    
    % Update variables
    iter = iter + 1;
    prev_norm = curr_norm;
  
    
end





% Calculation of the norm
norm_G=0.5*norm(sol,2)+gamma*norm_sumg(sol,param.G);

% Log after the calculous of the prox
if param.verbose >= 1
    fprintf('  prox_sumG: Sum_||G(x)|| = %e\n', norm_G);


    % Stopping criterion
    fprintf(' %i iterations\n', iter);
    fprintf(' Stopping criterion: %s \n\n', crit);

end


info.algo=mfilename;
info.iter=iter;
info.final_eval=norm_G;
info.crit=crit;
info.time=toc(t1);
end
