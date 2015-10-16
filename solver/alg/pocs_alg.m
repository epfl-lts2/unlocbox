function s = pocs_alg()
   s.name = 'POCS';
   s.initialize = @(x_0, fg, Fp, param) pocs_initialize(x_0,fg,param);
   s.algorithm = @(x_0, fg, Fp, sol, s, param) pocs_algorithm(Fp, sol, s);
   s.finalize = @(x_0, fg, Fp, sol, s, param) sol;

end

function [sol, s, param] = pocs_initialize(x_0,fg,param)
    


    sol = x_0;
    s = {};
    if fg.beta
        error('Beta = 0! This solver requires only projections functions.');
    end
end


function [sol, s] = pocs_algorithm(Fp, sol, s)
    
    for ii = 1 : length(Fp)
       sol = Fp{ii}.prox(sol,0);
    end

end
