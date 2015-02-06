function [sol, info,objective] = sdmm(F, param)
%SDMM Simultaneous-direction method of multipliers algorithm
%   Usage: sol = sdmm(F,param);
%          sol = sdmm(F);
%          [sol,info,objective] = sdmm(...);
%
%   Input parameters:
%         F     : Array of function to minimize
%         param : Optional parameter
%   Output parameters:
%         sol   : Solution
%         info : Structure summarizing informations at convergence
%         objective: vector (evaluation of the objectiv function each iteration)
%
%   `sdmm`, from simultaneous-direction method of multipliers solves:
% 
%   .. sol = argmin sum(f_i( L_i x))
%
%   .. math::  sol = \min_x \sum_i f_i( L_i x)
%
%   where $x$ belong to $R^N$, $L_i$ are linear operators and $x_i$ are the 
%   minimization variables.
%
%   * *F* is a cellarray of structure representing the functions.
%
%     In the function `F{i}`, there have to be:
%
%     * `F{i}.eval(x_i)` : an operator to evaluate the function
%     * `F{i}.prox(x_i, gamma)` : an operator to evaluate the prox of the function
%     * `F{i}.x0` : vector of initial value 
%
%     Optionally you can define
%
%     * `F{i}.L`  : linear operator, matrix or operator (default identity)
%     * `F{i}.Lt` : adjoint of linear operator, matrix or operator (default identity)
%
%
%   * *param* is a Matlab structure containing the following fields:
%
%     General parameters:
%
%     * *param.tol* : is stop criterion for the loop. The algorithm stops if
%
%       ..  max_i ||  y_i(t) - y_i(t-1) ||  / ||y(t) ||< tol,
%      
%       .. math:: \max_i \frac{ \| y_i(t) - y_i(t-1)\| }{ \|y_i(t)\|} < tol,
%
%       where  $y_i(t)$ are the dual variable of function *i* at itertion *t*
%       by default, `tol=10e-4`.
%
%       Warning! This stopping criterion is different from other solver!
%
%     * *param.maxit* : is the maximum number of iteration. By default, it is 200.
% 
%     * *param.verbose* : 0 no log, 1 print main steps, 2 print all steps.
%
%     * *param.gamma* : convergence parameter (default 1)
%
%     * *param.abs_tol* : If activated, this stopping critterion is the
%       objectiv function smaller than *param.tol*. By default 0.
%
%     * *param.Qinv* : Inverted Q matrix. $Qinv=Q^{-1}$ with:
%
%       .. Q = sum_i(L_i^t( L_i x))
%
%       .. math::  Q = \sum_i L_i^T ( L_i x)
%
%       By default, Qinv is the identity matrix divided by the number of
%       functions.
%
%       This parameter can be given in a matrix form or in a linear operator
%       form.
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
%   See also:  admm forward_backward douglas_rachford
%
%   Demos:  demo_sdmm
%
%   References: combettes2007douglas combettes2011proximal


% Author:  Nathanael Perraudin
% Date: fev 23 2012
%

% Start the time counter
t1 = tic;

% test the evaluate function
[F] = test_eval(F);

% number of function
    m = size(F,2);

% Optional input arguments
if nargin<2, param=struct; end

if ~isfield(param, 'tol'), param.tol=10e-4 ; end
if ~isfield(param, 'maxit'), param.maxit=200; end
if ~isfield(param, 'verbose'), param.verbose=1 ; end
if ~isfield(param, 'lambda'), param.lambda=1 ; end
if ~isfield(param, 'gamma'), param.gamma=1 ; end
if ~isfield(param, 'Qinv'), param.Qinv=@(x) x./m; end

if isa(param.Qinv,'numeric')
    param.QinvOp= @(x) param.Qinv*x;
else
	param.QinvOp= param.Qinv;
end

for i=1:m
    if ~isfield(F{i}, 'L'), F{i}.L=@(x) x; end
    if ~isfield(F{i}, 'Lt'), F{i}.Lt=@(x) x; end
    if size(F{i}.L,1)==0, F{i}.L=eye(length(F{i}.x0)) ; end
    if size(F{i}.x0,2)>size(F{i}.x0,1), F{i}.x0=F{i}.x0'; end
    F{i}.y_n=F{i}.x0;
    F{i}.z_n=F{i}.x0;
    F{i}.y_old=F{i}.x0;
    % Check how the matrix is given
    if isa(F{i}.L,'numeric')
        F{i}.OpL= @(x) F{i}.L*x;
    else
        F{i}.OpL= F{i}.L;
    end
    
    if isa(F{i}.Lt,'numeric')
        F{i}.OpLt= @(x) F{i}.Lt*x;
    else
        F{i}.OpLt= F{i}.Lt;
    end
end

% Initialization
x_n=F{1}.OpLt((F{1}.y_n-F{1}.z_n));
for i=2:m
    x_n=x_n+F{i}.OpLt((F{i}.y_n-F{i}.z_n));
end
x_n=param.QinvOp(x_n);

[~,~,prev_norm,iter,objective,~] = convergence_test(gen_norm(x_n,F));

% Main loop
while 1
    
    %
    if param.verbose >= 1
        fprintf('Iteration %i:\n', iter);
    end
    
    
    reldual=0;
    % Algorithm
    for i=1:m
        s_n=F{i}.OpL(x_n);
        F{i}.y_old=F{i}.y_n; % for convergence criterion
        F{i}.y_n=F{i}.prox(s_n+F{i}.z_n,param.gamma);
        F{i}.z_n=F{i}.z_n+s_n-F{i}.y_n ;% updates
        temp=norm(F{i}.y_old-F{i}.y_n)/norm(F{i}.y_n);
        if temp > reldual
            reldual=temp;
        end
    end
    x_n=F{1}.OpLt((F{1}.y_n-F{1}.z_n));
    for i=2:m
        x_n=x_n+F{i}.Lt((F{i}.y_n-F{i}.z_n));
    end
    x_n=param.QinvOp(x_n);
    sol=x_n; 

    
    % Global stopping criterion
    curr_norm = gen_norm(sol,F);  
    [~,rel_norm,prev_norm,iter,objective,crit] = convergence_test(...
            curr_norm,prev_norm,iter,objective,param);
    [sol,param] = post_process(sol, iter, curr_norm, prev_norm, objective, param);
    if reldual<param.tol
        break;
    end
    if iter>=param.maxit
        break;
    end
    if param.verbose >= 1
        fprintf([' ||f|| = %e, rel_norm = %e\n Maximum relative ' ...
            'distance of dual variable: %e\n'], ...
            curr_norm, rel_norm, reldual);
    end

    

    
end

% Log
if param.verbose>=1
    %Print norm
    fprintf('\n Solution found:\n');
    fprintf(' Final relative norm: %e\n', rel_norm );
    
    
    % Stopping criterion
    fprintf(' %i iterations\n', iter);
    fprintf(' Stopping criterion: %s \n\n', crit);
    
end

info.algo=mfilename;
info.iter=iter;
info.final_eval=curr_norm;
info.crit=crit;
info.time=toc(t1);
info.rel_norm=rel_norm;

end

function n=gen_norm(x_n,F)
% number of function
    m = size(F,2);
    n=0;
    for i=1:m
       n=n+F{i}.eval(F{i}.OpL(x_n)); 
    end
end
