function [stop, crit, s, iter, info] = convergence_test(sol, s, iter, fg, Fp, info, param)
%CONVERGENCE_TEST Test the convergence of an algorithm
%   Usage: stop = convergence_test(curr_norm,prev_norm,objective,param);
%
%   This function test the convergence of an algorithms and performs
%   some updates. Additionally, it handles verbosity during the iterations.
%
%   At convergence, the flag stop is set to one.
%

%

%   Nathanael Perraudin
%   Date: 16 oct 2015

global FS

%% Optional input arguments

% Param
if ~isfield(param, 'stop_box'), param.stop_box = 0; end

if isfield(param, 'use_dual') || isfield(param, 'abs_tol')
    warning('You are using the old method... Please update')
end

% Set algorithm name
if ~isfield(param, 'alg'), 
    txt = dbstack(); 
    param.alg = txt(2).name ; 
end

% Starting keyboard stop
if ~kbstop('lauched')
    kbstop('init');
end

% Use the stopbox
if param.stop_box
    if iter<=1
        if isstruct(FS)
            try  %#ok<TRYNC>
                FS.close();
            end
        end
        FS = stopstruct(param.alg,'Stop the algorithm') ;
    end
end


stop = 0;
crit= 'NOT_DEFINED';

switch lower(param.stopping_criterion)
    case 'rel_norm_obj'
        curr_eval = eval_objective(fg,Fp,sol,s,param);
        info.objective(iter+1) = curr_eval;
        rel_eval = abs(curr_eval - info.objective(iter))/(curr_eval + eps);
        info.rel_eval(iter) = rel_eval;
        if (abs(rel_eval) < param.tol) && iter>1
            crit = 'REL_NORM_OBJ';
            stop=1;   
        end

    case 'rel_norm_primal'
        rel_norm_primal = norm(sol(:)-s.prev_sol(:))/norm(sol(:));
        info.rel_norm_primal(iter) = rel_norm_primal;
        if (abs(rel_norm_primal) < param.tol) && iter>1
            crit = 'REL_NORM_PRIMAL';
            stop=1;   
        end
        s.prev_sol = sol;

    case 'rel_norm_dual'
        [rel_norm_dual, s.dual_var_old] = eval_dual_var(s.dual_var,s.dual_var_old);
        info.rel_norm_dual(iter) = rel_norm_dual;
        if (abs(rel_norm_dual) < param.tol) && iter>1
            crit = 'REL_NORM_DUAL';
            stop=1;   
        end   
    case 'rel_norm_primal_dual'
        rel_norm_primal = norm(sol(:)-s.prev_sol(:))/norm(sol(:));
        info.rel_norm_primal(iter) = rel_norm_primal;
        [rel_norm_dual, s.dual_var_old] = eval_dual_var(s.dual_var,s.dual_var_old);
        info.rel_norm_dual(iter) = rel_norm_dual;
        if (abs(rel_norm_primal) < param.tol) && ...
                (abs(rel_norm_dual) < param.tol) && iter>1
            crit = 'REL_NORM_PRIMAL_DUAL';
            stop=1;   
        end
        s.prev_sol = sol;
    case 'obj_increase'
        curr_eval = eval_objective(fg,Fp,sol,s,param);
        info.objective(iter+1) = curr_eval;
        if curr_eval >= info.objective(iter)
            crit = 'OBJ_INCREASE';
            stop = 1;   
        end
    case 'obj_threshold'
        curr_eval = eval_objective(fg,Fp,sol,s,param);
        info.objective(iter+1) = curr_eval;
        if (curr_eval < param.tol) && iter>1
            crit = 'OBJ_THRESHOLD';
            stop=1;    
        end        
    otherwise
        error('Unknown stopping criterion!')
end

%% Handling debug mode
% We compute what has not been computed yet
if param.debug_mode
    switch lower(param.stopping_criterion)
        case 'rel_norm_obj'
            rel_norm_primal = norm(sol(:)-s.prev_sol(:))/norm(sol(:));
            info.rel_norm_primal(iter) = rel_norm_primal;
            if isfield(s,'dual_var')
                [rel_norm_dual, s.dual_var_old] = ...
                    eval_dual_var(s.dual_var,s.dual_var_old);
                info.rel_norm_dual(iter) = rel_norm_dual;                
            end

        case 'rel_norm_primal'
            curr_eval = eval_objective(fg,Fp,sol,s,param);
            info.objective(iter+1) = curr_eval;
            info.rel_eval(iter) = ...
                abs((info.objective(iter+1)-info.objective(iter))/ ...
                info.objective(iter+1));

            if isfield(s,'dual_var')
                [rel_norm_dual, s.dual_var_old] = ...
                    eval_dual_var(s.dual_var,s.dual_var_old);
                info.rel_norm_dual(iter) = rel_norm_dual;                
            end

        case 'rel_norm_dual'
            curr_eval = eval_objective(fg,Fp,sol,s,param);
            info.objective(iter+1) = curr_eval;
            info.rel_eval(iter) = ...
                abs((info.objective(iter+1)-info.objective(iter))/ ...
                info.objective(iter+1));

            rel_norm_primal = norm(sol(:)-s.prev_sol(:))/norm(sol(:));
            info.rel_norm_primal(iter) = rel_norm_primal;
            s.prev_sol = sol;
            
        case 'rel_norm_primal_dual'
            curr_eval = eval_objective(fg,Fp,sol,s,param);
            info.objective(iter+1) = curr_eval;
            info.rel_eval(iter) = ...
                abs((info.objective(iter+1)-info.objective(iter))/ ...
                info.objective(iter+1));

        case 'obj_increase'
            rel_norm_primal = norm(sol(:)-s.prev_sol(:))/norm(sol(:));
            info.rel_norm_primal(iter) = rel_norm_primal;
            info.rel_eval(iter) = ...
                abs((info.objective(iter+1)-info.objective(iter))/ ...
                info.objective(iter+1));

            if isfield(s,'dual_var')
                [rel_norm_dual, s.dual_var_old] = ...
                    eval_dual_var(s.dual_var,s.dual_var_old);
                info.rel_norm_dual(iter) = rel_norm_dual;                
            end
        case 'obj_threshold'
            rel_norm_primal = norm(sol(:)-s.prev_sol(:))/norm(sol(:));
            info.rel_norm_primal(iter) = rel_norm_primal;
            info.rel_eval(iter) = ...
                abs((info.objective(iter+1)-info.objective(iter))/ ...
                info.objective(iter+1));
            if isfield(s,'dual_var')
                [rel_norm_dual, s.dual_var_old] = ...
                    eval_dual_var(s.dual_var,s.dual_var_old);
                info.rel_norm_dual(iter) = rel_norm_dual;                
            end     
        otherwise
            error('Unknown stopping criterion!')
    end    
end

    
%% Handle verbosity
if param.verbose >= 2
    if param.debug_mode
        fprintf('  f(x^*) = %e, rel_eval = %e\n', ...
            curr_eval, info.rel_eval(iter));       
        fprintf('   Relative norm of the primal variables: %e\n',...
                    rel_norm_primal);
        if isfield(s,'dual_var')
        	fprintf('   Relative norm of the dual variables: %e\n',...
            	rel_norm_dual);  
        end
    else
        switch lower(param.stopping_criterion)
            case 'rel_norm_obj'   
                fprintf('  f(x^*) = %e, rel_eval = %e\n', ...
                    curr_eval, rel_eval);                   
            case 'rel_norm_primal'
                fprintf('   Relative norm of the primal variables: %e\n',...
                    rel_norm_primal);
            case 'rel_norm_dual'
                fprintf('   Relative norm of the dual variables: %e\n',...
                    rel_norm_dual);                
            case 'rel_norm_primal_dual'
                fprintf('   Relative norm of the primal variables: %e\n',...
                    rel_norm_primal);
                fprintf('   Relative norm of the dual variables: %e\n',...
                    rel_norm_dual);  
            case 'obj_increase'
                fprintf('  f(x^*) = %e, prev_it:  %e\n', ...
                    curr_eval, info.objective(iter));   
            case 'obj_threshold'
                 fprintf('  f(x^*) = %e\n', ...
                    curr_eval);      
            otherwise
                error('Unknown stopping criterion!')
        end
    end
end

%% Stopping if too many iteration

if iter >= param.maxit
    crit= 'MAX_IT';
    stop=1;  
end

%% Stopping if the user press crtl +d
if  kbstop() || (param.stop_box && FS.stop())
    crit= 'USER';
    stop = 1;
end

if stop
    kbstop('stop');
    if param.stop_box

        try %#ok<TRYNC>
            FS.close() ;  % Clear up the box
        end

        try %#ok<TRYNC>
            clear FS ;    % this structure has no use anymore 
        end
    end
end


%% Performed updates
if ~stop
    iter=iter+1;
end




end


function curr_eval = eval_objective(fg,Fp,sol,s,param)
        curr_eval = eval_function(fg,Fp,sol,s,param);
        if ~(numel(curr_eval)==1)
            error('A least, one of your evaluation function does not return a scalar.')
        end

%         if isa(curr_norm,'gpuArray')
%             s.objective(iter)=gather(curr_norm);
%         else
%             s.objective(iter)=curr_norm;
%         end
end


function [rel_norm, x_old] =  eval_dual_var(x1,x2)

    if nargin<2
        if iscell(x1)
            x2 = cell(length(x1),1);
            for ii = 1:length(x1)
                x2{ii} = 0;
            end
        else
            x2 = 0;
        end
    end
    
    if iscell(x1)
        rel_norm = 0;
        for ii = 1:length(x1)
            tmp = norm(x1{ii}(:)-x2{ii}(:))/(norm(x1{ii}(:))+eps);
            if tmp>rel_norm
                rel_norm = tmp;
            end
        end
    else
        rel_norm = norm(x1(:)-x2(:))/(norm(x1(:))+eps);
    end


    x_old = x1;
end


