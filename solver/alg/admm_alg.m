function s = admm_alg()
   s.name = 'ADMM';
   s.initialize = @(x_0, fg, Fp, param) admm_initialize(x_0,fg,Fp,param);
   s.algorithm = @(x_0, fg, Fp, sol, s, param) admm_algorithm(Fp, s, param);
   s.finalize = @(x_0, fg, Fp, sol, s, param) sol;

end

function [sol, s, param] = admm_initialize(x_0,fg,Fp,param)
    
    if ~isfield(param, 'L'), param.L=@(x) x; end
    if isa(param.L,'numeric')
       s.OpL= @(x) param.L*x;
    else
       s.OpL= param.L;
    end
    s.x_n = {};
    s.x_n{2}{1} = s.OpL(x_0);
    s.u_n = zeros(size(s.x_n{2}{1}));
    s.dual_var = s.x_n{2}{1};
    sol = x_0;
    if fg.beta
        error('ADMM needs only function with proximal operators')
    end
    
    if ~(numel(Fp)==2)
        error('ADMM needs exactly 2 functions')
    end
    
    param.abs_tol = 1;
    param.use_dual = 1;
   
    
end


function [sol, s] = admm_algorithm(Fp, s, param)
    
    s.x_n{1} = Fp{1}.prox_ad(s.x_n{2}{1} - s.u_n,param.gamma);
    s_n= s .OpL(s.x_n{1}{1});
    s.x_n{2} = Fp{2}.prox_ad(s_n+s.u_n,param.gamma);
    s.dual_var = s.x_n{2}{1};

    
    s.u_n = s.u_n+s_n-s.x_n{2}{1} ;% updates
    sol = s.x_n{1}{1}; 
       
    
end
