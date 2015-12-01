function [sol, info] = solvep(x_0, F, param)
%SOLVEP solve a minimization problem
%   Usage: sol = solvep(x_0, F, param);
%          sol = solvep(x_0, F);
%          [sol,infos,objectiv] = solvep(...);
%
%   Input parameters:
%         x_0   : Starting point of the algorithm
%         F     : array of function to minimize (structure)
%         param : Optional parameter
%   Output parameters:
%         sol   : Solution
%         info  : Structure summarizing informations at convergence
%
%   `solvep` solves:
%
%   .. sol = argmin sum_i fi(x)     for x belong to R^N
%
%   .. math::  sol = arg \min_x |sum_i f_i(x)  \hspace{1cm} for \hspace{1cm}  x\in R^N
%
%   where *x* is the variable.
%
%   *x_0* is the starting point of the algorithm. A good starting point
%   could significantly reduce the computation time
%
%   *F* is an array of structure representing convex function to be
%   minimized. Those Functions
%   can be minimized thanks to their gradient (if they are differentiable)
%   or thanks to their proximal operator. As a result the algorithm will
%   needs some of those operators. The easiest way to define a 
%   function *f1* is to create a struture with the fields f1.eval, f1.grad
%   and f1.prox. Those field all contatains an inline function that
%   compute respectively the evaluation of the function itself, the
%   gradient or the proximal operator. Depending on the solver, not all
%   this operators are necessary. See each solver documentation for
%   details. When three functions are defined, F = {f1, f2, f3}.
%   
%   *param* a Matlab structure containing the following fields:
%
%   * *param.gamma* : is the step size. Watch out, this parameter is
%     bounded. It should be below $1/\beta$ (*f2* is $\beta$ Lipchitz
%     continuous). By default, it is computed with the lipschitz constant
%     of all smooth functions.
%
%   * *param.tol* : Tolerance to stop iterating. Please see
%     *param.stopping_criterion*. (Default 1e-4).
%
%   * *param.algo* : solver used for the problem. Determined
%     automatically with the functions in *f*.
%
%   * *param.stopping_criterion* : is stopping criterion to end the
%     algorithm. Possible values are:
%
%     * 'rel_norm_obj' : Relative norm of the objective function.
%     * 'rel_norm_primal' : Relative norm of the primal variables. 
%     * 'rel_norm_dual' : Relative norm of the dual variables. 
%     * 'rel_norm_primal_dual' : Relative norm of the primal and the dual
%       variables.  
%     * 'obj_increase' : Stops when the objective function starts
%       increasing or stay equal.
%     * 'obj_threshold' : Stops when the objective function is below a
%       threshold. The threshold is set in *param.tol*.
%
%     For the 'rel_norm' stopping criterion, the algorithm end if 
%
%       ..  ||  n(t) - n(t-1) ||_2  /  || n(t) ||_2 < tol,
%      
%       .. math:: \frac{ \| n(t) - n(t-1) \|_2 }{\| n(t)\|_2} < tol,
%
%       where  $n(t)$ is the objective function, the primal or the dual
%       variable at iteration *t*. 
%
%   * *param.maxit* : is the maximum number of iteration. By default, it is
%     200. 
% 
%   * *param.verbose* : 0 no log, 1 print main steps, 2 print all steps.
%
%   * *param.debug_mode* : Compute all internal convergence parameters.
%     Activate this option for debugging
%
%   *info* is a Matlab structure containing the following fields:
%
%   * *info.algo* : Algorithm used
%
%   * *info.iter* : Number of iteration
%
%   * *info.time* : Time of exectution of the function in sec.
%
%   * *info.crit* : Stopping critterion used 
%
%   Additionally, depending on the stopping critterion, the structure
%   *info* also contains:
%
%   * *info.objective* : Value of the objective function
%
%   * *info.rel_norm_obj* : Relative norm of the objective function.
%
%   * *info.rel_norm_primal* : Relative norm of the primal variable.
%
%   If the flag *param.debug_mode* is activated, the previous quantity are
%   always computed. Moreover, for solver using dual variable, *info* also
%   contains:
%
%   * *info.rel_norm_dual* : Relative norm of the dual variable.
%
%   * *info.dual_var* : Final dual variables.
%   


% Author: Nathanael Perraudin
% Date: 16 oct 2015
% Testing: test_solver

% Start the time counter
t1 = tic;

% Optional input arguments
if nargin<3, param=struct; end

if nargin<2
    error('Not enough input arguments')
end

if nargout>3
    error('The third argument objective has been moved in info.objective')
end


if ~isfield(param, 'tol'), param.tol=1e-4 ; end
if ~isfield(param, 'maxit'), param.maxit=200; end
if ~isfield(param, 'verbose'), param.verbose=1 ; end
if ~isfield(param, 'lambda'), param.lambda=0.99 ; end
if ~isfield(param, 'fast_eval'), param.fast_eval = 0  ; end
if ~isfield(param, 'debug_mode'), param.debug_mode = 0 ; end



% test the input for eval
if ~iscell(F)
    F = {F};
end
F = test_eval(F);


% test input for grad and prox
[Fg,Fp] = prepare_function(F,param);
if ~isfield(param, 'gamma'), 
    if numel(Fg)
        param.gamma = compute_gamma(Fg);
    else
        param.gamma = 1; 
    end
else
    if param.verbose >= 1
        fprintf('The time step is set manually to : %g\n', param.gamma);
    end    
end

% Select the algorithm
if ~isfield(param, 'algo'), param.algo = select_solver(Fg,Fp)  ; end

% Select the stopping critterion
if ~isfield(param, 'stopping_criterion')
    param.stopping_criterion = select_stopping_criterion(param.algo);
end;

algo = get_algo(param.algo);

% Transform all smooth functions into one function.
fg = add_smooth_functions(Fg);

if param.verbose>1, 
    fprintf(['Algorithm selected: ', algo.name,' \n']);
end

[sol,s,param] = algo.initialize(x_0, fg, Fp, param);


% Initialization
[info, iter, s] = initialize_convergence_variable(sol, s, fg, Fp, param);

% Main loop
while 1

    if param.verbose >1
        fprintf('Iteration %i:\n', iter);
    end
    
    [sol, s] = algo.algorithm(x_0, fg, Fp, sol, s, param);
    
    [stop, crit, s, iter, info] = ...
        convergence_test(sol, s, iter, fg, Fp, info, param);
    [sol, param] = post_process(sol, iter, info, param);
    if stop
        break;
    end
    
end

sol = algo.finalize(x_0, fg, Fp, sol, s, param);

summary_print(s,info,iter,algo,crit,param);


info.algo = algo.name;
info.iter = iter;
info.crit = crit;
info.time = toc(t1);

% Return the dual variable if availlable
if isfield(s,'dual_var')
    info.dual_var = s.dual_var;
end

end





function gamma = compute_gamma(Fg)
    beta = 0;
    for ii = 1:length(Fg)
        beta = beta + Fg{ii}.beta;
    end
    if beta == 0;
        gamma = 1;
    elseif beta >0
        gamma = 1/beta;
    else
        error('Error in the libschitz constant!')
    end
end

function algo = get_algo(name)
    if isstruct(name)
        algo = name;
    elseif ischar(name)
        algo = algoname(name);
    else
        error('The algorithm is not a string and not a struct!')
    end
end
