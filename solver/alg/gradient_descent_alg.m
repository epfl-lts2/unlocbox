function s = gradient_descent_alg()
   s.name = 'GRADIENT_DESCENT';
   s.initialize = @(x_0, fg, Fp, param) gradient_descent_initialize(x_0,fg,Fp,param);
   s.algorithm = @(x_0, fg, Fp, sol, s, param) gradient_descent_algorithm(fg, sol, s, param);
   s.finalize = @(x_0, fg, Fp, sol, s, param) sol;

end

function [sol, s, param] = gradient_descent_initialize(x_0,fg,Fp,param)
    

%     s.u_n = x_0;
%     s.tn = 1;
%     sol = s.u_n - param.gamma*fg.grad(s.u_n);
    sol = x_0;
    s = {};
    
    if numel(Fp)>0
        error('This solver can not be used to optimize non smooth functions')
    end
    
    if ~fg.beta
        error('Beta = 0! This solver requires a smooth term.');
    end
end


function [sol, s] = gradient_descent_algorithm(fg, sol, s, param)
    
    
%     tn1 = (1 + sqrt(1+4*s.tn^2))/2;
%     tmp = sol + (s.tn-1)/tn1*(s.u_n - sol);
%     s.u_n = sol;
%     sol = tmp - param.gamma*fg.grad(tmp);
%     s.tn = tn1;
%     
    sol = sol-param.gamma*fg.grad(sol);


end
