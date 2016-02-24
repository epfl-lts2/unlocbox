function [sol,iter,objectiv] = ppxa_old(x_0,F,param)
%PPXA_OLD Parallel Proximal algorithm
%   Usage: sol = ppxa_old(x_0,F,param);
%          sol = ppxa_old(x_0,F);
%          [sol,iter] = ppxa_old(...);
%
%   Input parameters:
%         x_0   : Starting point of the algorithm
%         F     : Array of function to minimize
%         param : Optional parameter
%   Output parameters:
%         sol   : Solution
%         iter  : Current iteration
%
%   `ppxa`, derived from the Douglas-Rachford algorithm, solves
% 
%   .. sol = argmin sum(W_i*f_i(x))
%
%   .. math::  sol = \min_x \sum_i W_i f_i(x)
%
%   for *x* in $R^N$, where *x* is the variable and *x_0* is the starting point.
%
%   * *F* is an array of structures representing functions containing 
%     operators inside and eventually the norm. If no norm is defined, the 
%     $l^1$ norm is used. The prox: *F(i).prox* and norm: *F(i).eval* are defined 
%     in the same way as in the Forward-backward and Douglas-Rachford algorithms
%
%   * *param* a Matlab structure containing the following fields:
%
%     General parameters:
%
%     * *param.W* : the weight (all equal by default)
%
%     * *param.tol* : is stop criterion for the loop. The algorithm stops if
%
%       ..  (  n(t) - n(t-1) )  / n(t) < tol,
%      
%       .. math:: \frac{  n(t) - n(t-1) }{ n(t)} < tol,
%
%       where  $n(t) = \sum W_i*f_i(x)$ is the objective function at iteration *t*
%       by default, `tol=10e-4`.
%
%     * *param.maxit* : is the maximum number of iteration. By default, it is 200.
% 
%     * *param.verbose* : 0 no log, 1 print main steps, 2 print all steps.
%
%     * *param.lambda* : is the weight of the update term.
%
%     * *param.gamma* : convergence parameter (default 1)
%
%     * *param.abs_tol* : If activated, this stopping critterion is the
%       objectiv function smaller than *param.tol*. By default 0.
%
%   See also:  sdmm, admm, generalized_forward_backward
%
%   Demos:  demo_ppxa
%
%   References:  combettes2011proximal

% Author:  Nathanael Perraudin
% Date: Oct 19 2012



% number of function
m = size(F,2);

    % Optional input arguments
    if nargin<3, param=struct; end

    if ~isfield(param, 'gamma'), param.gamma = 1; end
    if ~isfield(param, 'tol'), param.tol=10e-4 ; end
    if ~isfield(param, 'maxit'), param.maxit=200; end
    if ~isfield(param, 'verbose'), param.verbose=1 ; end
    if ~isfield(param, 'lambda'), param.lambda=0.99 ; end
    if ~isfield(param, 'W'), param.W=ones(m,1)/m ; end
    
    W=param.W;
    
    
   % Reshape x if vector line
    if (size(W,1)>size(W,2))
        W=W';
    end
    
    test_gamma(param.gamma);
    test_sum(W);
    
    
    
    sz = size(x_0,1);
    w = repmat(W',1,sz);
    y = zeros(m,sz);
    p = zeros(m,sz);
    
    
    %Initilisation
    for i = 1:m
         y(i,:) = x_0;
    end

    x = sum(w.*y);

    
    % outerloop
    curr_norm = 0;
    for i = 1:m
       curr_norm = F(i).eval(x)+curr_norm;
    end
    [~,~,prev_norm,iter,objectiv,~] = convergence_test_old(curr_norm);
    
    while 1
        
        if param.verbose >= 1
            fprintf('Iteration %i:\n', iter);
        end
        
        % parallel proximal
        % compute updated prox
        for i = 1:m
            temp = F(i).prox(y(i,:)', param.gamma);
            p(i,:) = temp(:);
        end
        pn = sum(w.*p);

        param.lambda = update_lambda(param.lambda,iter);

        % update y
        y  = y + param.lambda*(2*repmat(pn,m,1)-repmat(x,m,1)-p);

        % update x
        x = x + param.lambda * (pn - x);
        
        % update solution & relative norm
        sol = x';
        curr_norm = 0;
        for i = 1:m
            curr_norm = F(i).eval(sol)+curr_norm;
        end
                
        [stop,rel_norm,prev_norm,iter,objectiv,crit] = convergence_test_old(curr_norm,prev_norm,iter,objectiv,param);
        if stop
            break;
        end
        if param.verbose >= 1
            fprintf('Current objectiv function : %e \t relative norm : %e \n', curr_norm, rel_norm);
        end
        
    end
    
    if param.verbose>=1
        % Stopping criterion
        fprintf(' %i iterations\n', iter);
        fprintf(' Stopping criterion: %s \n\n', crit);
    end
    
end

function lambda = update_lambda(lambda,n)
    lambda = lambda;
end


function test_gamma(gamma)
    if gamma <= 0 
        fprintf('Warning : gamma is not > 0\n');
    end
end


function test_sum(W)
    if sum(W) == 1
    else
        fprintf('Warning : sum W is not equal to 1');
    end
end
