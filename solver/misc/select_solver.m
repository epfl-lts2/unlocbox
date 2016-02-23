function solver = select_solver(Fg,Fp)
%SELECT_SOLVER This function choose a default solver
%   Usage: solver = select_solver(Fg,Fp)
%
%   This function choose a default solver depending on the function Fg and
%   Fp. Fg is cell array of smooth functions and Fp a cell array of non
%   smooth functions.

    n = numberofL(Fp);
    
    if numel(Fg) % There are smooth functions
        if n>2
            error('Sorry, no solver is able to solve your problem yet!')
        end
        if numel(Fp)==0
            solver = 'GRADIENT_DESCENT';
        elseif numel(Fp)==2 
            solver = 'FB_BASED_PRIMAL_DUAL';
        elseif (numel(Fp)<=2) && n
            solver = 'FB_BASED_PRIMAL_DUAL';            
        elseif numel(Fp)==1
            solver = 'FORWARD_BACKWARD';
        else
            solver = 'GENERALIZED_FORWARD_BACKWARD';
        end
    else % All function are non smooth
     
        if (numel(Fp)<=2) && (n == 1)
            solver = 'CHAMBOLLE_POCK';    
        elseif (n>1)
            solver = 'SDMM';
        elseif numel(Fp) == 2
            solver = 'DOUGLAS_RACHFORD';
        else
            solver = 'PPXA';
        end
    end
        
end




function n = numberofL(Fp)
% Return the number of functions with a linear opeartor inside
n = 0;
    for ii = 1:length(Fp)
        if isfield(Fp{ii},'L') 
            n = n + 1;
            if ~isfield(Fp{ii},'Lt')
                warning('You did not define the Lt operator!');
            end
        end
    end
end