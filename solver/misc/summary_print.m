function [] = summary_print(s,info,iter,algo,crit,param)
%SUMMARY_PRINT Print the summary at the end of the optimization
%
%   This function print a summary at the end of the optimization depending
%   on the verbosity, the stopping critertion and the computed
%   statistiques.

%   Nathanael Perraudin
%   Date: 16 oct 2015

%% Handle verbosity
if param.verbose >= 2
    fprintf(['\n ',algo.name,':\n']);

    if param.debug_mode
        fprintf('  f(x^*) = %e, rel_eval = %e\n', ...
            info.objective(iter+1), info.rel_eval(iter));       
        fprintf('   Relative norm of the primal variables: %e\n',...
                    info.rel_norm_primal(iter));
        if isfield(s,'dual_var')
        	fprintf('   Relative norm of the dual variables: %e\n',...
            	info.rel_norm_dual(iter));  
        end
    else
        switch lower(param.stopping_criterion)
            case 'rel_norm_obj'   
                fprintf('  f(x^*) = %e, rel_eval = %e\n', ...
                    info.objective(iter+1), info.rel_eval(iter));                   
            case 'rel_norm_primal'
                fprintf('   Relative norm of the primal variables: %e\n',...
                    info.rel_norm_primal(iter));
            case 'rel_norm_dual'
                fprintf('   Relative norm of the dual variables: %e\n',...
                    info.rel_norm_dual(iter));                
            case 'rel_norm_primal_dual'
                fprintf('   Relative norm of the primal variables: %e\n',...
                    info.rel_norm_primal(iter));
                fprintf('   Relative norm of the dual variables: %e\n',...
                    info.rel_norm_dual(iter));  
            case 'obj_increase'
                fprintf('  f(x^*) = %e, prev_it:  %e\n', ...
                    info.objective(iter), info.objective(iter-1));   
            case 'obj_threshold'
                 fprintf('  f(x^*) = %e\n', ...
                    info.objective(iter));      
            otherwise
                error('Unknown stopping criterion!')
        end
    end
    % Stopping criterion
    fprintf(' %i iterations\n', iter);
    fprintf(' Stopping criterion: %s \n\n', crit);
elseif param.verbose>=1
%     if param.use_dual && isfield(s,'dual_var')
%         fprintf([algo.name,': Final dual relative evaluation = %e, it = %i, %s\n'], ...
%                         curr_eval, iter,crit);
%     else
%         fprintf([algo.name,': f(x^*) = %e, rel_eval = %e, it = %i, %s\n'], ...
%                         curr_eval, rel_eval, iter,crit);
%     end
        switch lower(param.stopping_criterion)
            case 'rel_norm_obj'   
                fprintf([algo.name,'  f(x^*) = %e, rel_eval = %e, it = %i, %s\n'], ...
                    info.objective(iter+1), info.rel_eval(iter), iter, crit);                   
            case 'rel_norm_primal'
                fprintf([algo.name,'   Relative norm of the primal variables: %e, it = %i, %s\n'],...
                    info.rel_norm_primal(iter), iter, crit);
            case 'rel_norm_dual'
                fprintf([algo.name,'   Relative norm of the dual variables: %e, it = %i, %s\n'],...
                    info.rel_norm_dual(iter), iter, crit);             
            case 'rel_norm_primal_dual'
                fprintf([algo.name,'   Rel primal: %e, rel dual %e, it = %i, %s\n'],...
                    info.rel_norm_primal(iter), info.rel_norm_dual(iter), iter, crit);
            case 'obj_increase'
                fprintf([algo.name,'  f(x^*) = %e, prev_it:  %e, it = %i, %s\n'], ...
                    info.objective(iter), info.objective(iter+1), iter, crit);  
            case 'obj_threshold'
                 fprintf([algo.name,'  f(x^*) = %e, it = %i, %s\n'], ...
                    info.objective(iter), iter, crit);  
            otherwise
                error('Unknown stopping criterion!')
        end
end
    
end







