function s = demo_forward_backward_alg()
%DEMO_FORWARD_BACKWARD_ALG Demonstration to define a personal solver
%   Usage : param.algo = demo_forward_backward_alg();
%
%   This function returns a structure containing the algorithm. You can
%   lauch your personal algorithm with the following::
%
%           param.algo = demo_forward_backward_alg();
%           sol = solvep(x0, {f1, f2}, param);
%

    % This function returns a structure with 4 fields:
    % 1) The name of the solver. This is used to select the solvers.
    s.name = 'DEMO_FORWARD_BACKWARD';
    % 2) A method to initialize the solver (called at the beginning)
    s.initialize = @(x_0, fg, Fp, param) ...
      forward_backward_initialize(x_0,fg,Fp,param);
    % 3) The algorithm itself (called at each iterations)
    s.algorithm = @(x_0, fg, Fp, sol, s, param) ...
      forward_backward_algorithm(fg, Fp, sol, s, param);
    % 4) Post process method (called at the end)
    s.finalize = @(x_0, fg, Fp, sol, s, param) sol;
    % The variables here are
    %   x_0 : The starting point
    %   fg  : A single smooth function 
    %         (if fg.beta == 0, no smooth function is specified) 
    %   Fp  : The non smooth functions (a cell array of structure)
    %   param: The structure of optional parameter
    %   s   : Intern variables or the algorithm
    %   sol : Current solution
    % Note that in s some field are protected: 
    %   * s.dual_var : cell array or array for the dual variable
    %   * s.dual_var_old : to store the variable of the previous iteration
    %   * s.prev_sol : previous solution
    % The field *s.dual_var* should be updated in this function. The other
    % should not be modified.
end

function [sol, s, param] = forward_backward_initialize(x_0,fg,Fp,param)

    % Handle optional parameter. Here we need a variable lambda.
    if ~isfield(param, 'lambda'), param.lambda=1 ; end

    % All intern variables are stored into the structure s
    s = struct;
    % *sol* is set to the initial points
    sol = x_0;
    
    if numel(Fp)>1
        error(['This solver can not be used to optimize',...
          ' more than one non smooth function']);
    end
    
    if ~fg.beta
        error('Beta = 0! This solver requires a smooth term.');
    end
    
end

function [sol, s] = forward_backward_algorithm(fg, Fp, sol, s, param)
    % The forward backward algorithm is done in two steps
    %  1) x_n = prox_{f, gamma} ( sol - gamma grad_fg(sol) )
    s.x_n = Fp{1}.prox( sol - param.gamma*fg.grad(sol), param.gamma);
    %  2) Updates
    %     sol = sol + lambda * (x_n -sol)
    sol = sol + param.lambda * (s.x_n - sol);
end