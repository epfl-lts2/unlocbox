function [info, iter, s] = initialize_convergence_variable(sol, s, fg, Fp, param)


% Test if the the dual variables are present
if strcmpi(param.test_type, 'rel_norm_primal_dual') ...
    || strcmpi(param.test_type, 'rel_norm_dual');
    if ~isfield(s,'dual_var')
        error('The dual variable are not ignitialized. Put them in s.dual_var');
    end
end


% In debut mode we inialize everything
if param.debug_mode
    info.objective = nan(param.maxit+1, 1);
    info.ogjective(1) = eval_function(fg, Fp, sol, s, param);
    s.prev_sol = sol;
    info.rel_norm_primal = nan(param.maxit, 1);
    info.rel_norm_dual = nan(param.maxit, 1);
    info.rel_eval = nan(param.maxit,1);
    if isfield(s,'dual_var')
    	s.dual_var_old = s.dual_var;
    end
else

    switch lower(param.test_type)
        case 'rel_norm_obj'
            info.objective = nan(param.maxit+1,1);
            info.rel_eval = nan(param.maxit,1);
            info.ogjective(1) = eval_function(fg, Fp, sol, s, param);
        case 'rel_norm_primal'
            s.prev_sol = sol;
            info.rel_norm_primal = nan(param.maxit, 1);
        case 'rel_norm_dual'
            s.dual_var_old = s.dual_var;
            info.rel_norm_dual = nan(param.maxit, 1);
        case 'rel_norm_primal_dual'
            s.dual_var_old = s.dual_var;
            s.prev_sol = sol;
            info.rel_norm_primal = nan(param.maxit, 1);
            info.rel_norm_dual = nan(param.maxit, 1);           
        case 'obj_increase'
            info.objective = nan(param.maxit+1,1);
            info.ogjective(1) = eval_function(fg, Fp, sol, s, param);
        case 'obj_threshold'
            info.objective = nan(param.maxit+1,1);
            info.ogjective(1) = eval_function(fg, Fp, sol, s, param);   
        otherwise
            error('Unknown stopping criterion!')
    end
end   


% We start at iteration 1
iter = 1;


end