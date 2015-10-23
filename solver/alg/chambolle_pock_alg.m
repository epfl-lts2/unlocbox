function s = chambolle_pock_alg()
   s.name = 'CHAMBOLLE_POCK';
   s.initialize = @(x_0, fg, Fp, param) chambolle_pock_initialize(x_0,fg,Fp,param);
   s.algorithm = @(x_0, fg, Fp, sol, s, param) chambolle_pock_algorithm(Fp, sol,s, param);
   s.finalize = @(x_0, fg, Fp, sol, s, param) sol;

end

function [sol, s, param] = chambolle_pock_initialize(x_0,fg,Fp,param)
    
    error('This solver contains a bug and is not working yet!')
    
    if ~isfield(param, 'tau'), param.tau=1 ; end
    if ~isfield(param, 'rho'), param.rho=1 ; end
    if ~isfield(param, 'L'), param.L=@(x) x; end
    if ~isfield(param, 'Lt'), param.Lt=@(x) x; end

    s.tau = param.tau;
    s.rho = param.rho;
    if isa(param.L,'numeric')
       s.OpL= @(x) param.L*x;
    else
       s.OpL= param.L;
    end

    if isa(param.Lt,'numeric')
       s.OpLt= @(x) param.Lt*x;
    else
       s.OpLt= param.Lt;
    end



    s.dual_var = s.OpL(x_0);
    sol = x_0;
    s.x_n = x_0;


    if fg.beta
        error('CHAMBOLLE POCK needs only function with proximal operators')
    end
    
    if ~(numel(Fp)==2)
        error('CHAMBOLLE POCK needs exactly 2 functions')
    end
    
   
    
end


function [sol, s] = chambolle_pock_algorithm(Fp, sol, s, param)
    
    % Algorithm
    s.dual_var = prox_adjoint( s.dual_var + s.rho *s.OpL(sol), s.rho, Fp{2});
    x_n_old = s.x_n;
    s.x_n = Fp{1}.prox( s.x_n + s.tau * s.OpLt(s.dual_var), s.tau);
    sol = s.x_n + param.gamma*(s.x_n - x_n_old);
  
end
