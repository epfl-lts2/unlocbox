function s = admm_alg()
   s.name = 'ADMM';
   s.initialize = @(x_0, fg, Fp, param) admm_initialize(x_0,fg,Fp,param);
   s.algorithm = @(x_0, fg, Fp, sol, s, param) admm_algorithm(Fp, s);
   s.finalize = @(x_0, fg, Fp, sol, s, param) sol;

end

function [sol, s, param] = admm_initialize(x_0,fg,Fp,param)
    
    if fg.beta
        error('ADMM needs only function with proximal operators')
    end
    
    if ~(numel(Fp)==2)
        error('ADMM needs exactly 2 functions')
    end
    s.x_n = {};
    
    if isfield(Fp{1},'proxL')
        s.ind = [1,2];
        L = Fp{2}.L;
    elseif isfield(Fp{2},'proxL')
        s.ind = [2,1];
        L = Fp{1}.L;        
    else
        L =@(x) x;
        s.ind = 0;
        if param.verbose
            warning('Warning! No operator L for admm!')
        end
    end
 
    if isa(L,'numeric')
       s.OpL= @(x) L*x;
    else
       s.OpL= L;
    end
%     s.x_n{2}{1} = s.OpL(x_0);
%     s.u_n = zeros(size(s.x_n{2}{1}));
%     s.dual_var = s.x_n{2}{1};    
    s.x_n{2} = s.OpL(x_0);
    s.u_n = zeros(size(s.x_n{2}));
    s.dual_var = s.x_n{2};
    sol = x_0;
    s.gamma = 1/param.gamma;

       
    
end


function [sol, s] = admm_algorithm(Fp, s)
    if s.ind(1)
%         s.x_n{1} = Fp{s.ind(1)}.prox_adL(s.x_n{2}{1} - s.u_n,s.gamma);
%         s_n= s.OpL(s.x_n{1}{1});
%         s.x_n{2} = Fp{s.ind(2)}.prox_ad(s_n+s.u_n,s.gamma);
        s.x_n{1} = Fp{s.ind(1)}.proxL(s.x_n{2} - s.u_n,s.gamma);
        s_n= s.OpL(s.x_n{1});
        s.x_n{2} = Fp{s.ind(2)}.prox(s_n+s.u_n,s.gamma);

    else
%         s.x_n{1} = Fp{1}.prox_ad(s.x_n{2}{1} - s.u_n,s.gamma);
%         s_n= s.OpL(s.x_n{1}{1});
%         s.x_n{2} = Fp{2}.prox_ad(s_n+s.u_n,s.gamma);        
        s.x_n{1} = Fp{1}.prox(s.x_n{2} - s.u_n,s.gamma);
        s_n= s.OpL(s.x_n{1});
        s.x_n{2} = Fp{2}.prox(s_n+s.u_n,s.gamma);        
    end
%     s.dual_var = s.x_n{2}{1};  
%     s.u_n = s.u_n+s_n-s.x_n{2}{1} ;% updates
%     sol = s.x_n{1}{1};     
    
    s.dual_var = s.x_n{2};
    s.u_n = s.u_n+s_n-s.x_n{2} ;% updates
    sol = s.x_n{1}; 
       
    
end
