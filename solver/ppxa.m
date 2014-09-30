function [sol, info, objective] = ppxa(x_0, F, param)
%PPXA Parallel Proximal algorithm
%   Usage: sol = ppxa(x_0, F, param);
%          sol = ppxa(x_0, F);
%          [sol, infos, objective] = ppxa(...); 
%
%   Input parameters:
%         x_0   : Starting point of the algorithm
%         F     : Array of function to minimize
%         param : Optional parameter
%   Output parameters:
%         sol   : Solution
%         info  : Structure summarizing informations at convergence
%         objective: vector (evaluation of the objective function each iteration)
%
%   `ppxa`, derived from the Douglas-Rachford algorithm, solves
% 
%   .. sol = argmin sum(W_i*f_i(x))
%
%   .. math::  sol = \min_x \sum_i W_i f_i(x)
%
%   for *x* in $R^N$, where *x* is the variable and *x_0* is the starting point.
%
%   * *F* is a cellarray of structures representing functions containing 
%     operators inside and eventually the norm. The prox: *F{i}.prox* and 
%     the function: *F{i}.eval* are defined in the same way as in the 
%     Forward-backward and Douglas-Rachford algorithms
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
%       objective function smaller than *param.tol*. By default 0.
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
%   * *info.final_eval* : Final evaluation of the objectivs functions
%
%   * *info.crit* : Stopping critterion used 
%
%   * *info.rel_norm* : Relative norm at convergence 
%
%
%   See also:  sdmm, admm, generalized_forward_backward
%
%   Demos:  demo_ppxa
%
%   References:  combettes2011proximal

% Author:  Nathanael Perraudin
% Date: Oct 19 2012

% Start the time counter
t1 = tic;

% test the evaluate function
[F] = test_eval(F);

% number of functions
m = size(F, 2);

    % Optional input arguments
    if nargin < 3, param = struct; end

    if ~isfield(param, 'gamma'), param.gamma = 1; end
    if ~isfield(param, 'tol'), param.tol = 10e-4; end
    if ~isfield(param, 'abs_tol'), param.abs_tol = 0; end
    if ~isfield(param, 'maxit'), param.maxit = 200; end
    if ~isfield(param, 'verbose'), param.verbose = 1; end
    if ~isfield(param, 'lambda'), param.lambda = 0.99; end
    if ~isfield(param, 'W'), param.W = ones(m, 1) / m; end

    
    W=param.W;
    
    
   % Reshape x if vector line
    if (size(W, 2) > size(W, 1))
        W = W';
    end
    
    test_gamma_ppxa(param.gamma);
    test_sum(W);
    
    % Create a table of scructure containing data
    data = [];
    for ii = 1:m
        data(ii).y = x_0;     %#ok<AGROW>
        data(ii).p = zeros(size(x_0));     %#ok<AGROW>
    end
    

    x = w_sum(W, data, 'y');

    
    % outerloop
    curr_norm = 0;
    for ii = 1:m
       curr_norm = F{ii}.eval(x) + curr_norm;
    end
    [~, ~, prev_norm, iter, objective, ~] = convergence_test(curr_norm);
    
    while 1
        
        if param.verbose >= 2
            fprintf('Iteration %i:\n', iter);
        end
        
        % parallel proximal
        % compute updated prox
        for ii = 1:m
            data(ii).p = F{ii}.prox(data(ii).y, param.gamma);%#ok<AGROW>
        end
        pn = w_sum(W,data,'p');



        % update y
        
        for ii = 1:m
            data(ii).y  = data(ii).y + param.lambda * (2*pn - x - data(ii).p); %#ok<AGROW>
        end

        % update x
        x = x + param.lambda * (pn - x);
        
        % update solution & relative norm
        sol = x;
        curr_norm = 0;
        for ii = 1:m
            curr_norm = F{ii}.eval(sol) + curr_norm;
        end
                
        [stop, rel_norm, prev_norm, iter, objective, crit] = convergence_test(curr_norm, prev_norm, iter, objective, param);
        [x, param] = post_process(sol, iter, curr_norm, prev_norm, param);
        if stop
            break;
        end
        if param.verbose >= 2
            fprintf('Current objective function : %e \t relative norm : %e \n', curr_norm, rel_norm);
        end
        
    end
    
    if param.verbose>=2
        % Print norm
        fprintf('\n Solution found:\n');
        if param.abs_tol
            fprintf(' Final norm: %e\n', curr_norm );
        else
            fprintf(' Final relative norm: %e\n', rel_norm );
        end
    
    
        % Stopping criterion
        fprintf(' %i iterations\n', iter);
        fprintf(' Stopping criterion: %s \n\n', crit);
    elseif param.verbose>=1
        fprintf('  Solution found: ||f|| = %e, rel_norm = %e, %s\n', ...
                    curr_norm, rel_norm,crit);
    end
    
info.algo = mfilename;
info.iter = iter;
info.final_eval = curr_norm;
info.crit = crit;
info.time = toc(t1);
info.rel_norm = rel_norm;

end



function test_gamma_ppxa(gamma)
    if gamma <= 0 
        fprintf('Warning : gamma is not > 0\n');
    end
end


function test_sum(W)
    if (sum(W) > 1+eps) || (sum(W) < 1-eps)
        fprintf('Warning : sum W is not equal to 1');
    end
end

function s = w_sum(w, data, l)
    s = zeros(size(data(1)));
    
    if l == 'p'
        for ii = 1:length(data)
            s = s + w(ii) * data(ii).p;
        end
    elseif l == 'y'
        for ii = 1:length(data)
            s = s + w(ii) * data(ii).y;
        end
    else
        error('Fatal error! Unknown field in data!')
    end
end