function s = douglas_rachford_alg()
   s.name = 'DOUGLAS_RACHFORD';
   s.initialize = @(x_0, fg, Fp, param) douglas_rachford_initialize(x_0,fg,Fp,param);
   s.algorithm = @(x_0, fg, Fp, sol, s, param) douglas_rachford_algorithm(Fp, sol, s, param);
   s.finalize = @(x_0, fg, Fp, sol, s, param) sol;

end

function [sol, s, param] = douglas_rachford_initialize(x_0,fg,Fp,param)

    if ~isfield(param, 'lambda'), param.lambda=1 ; end
    s.lambda = param.lambda;
    s.x_n = {};
    s.u_n = x_0;
    sol = x_0;
    if fg.beta
        error('Douglas rachford needs only function with proximal operators')
    end
    
    if ~(numel(Fp)==2)
        error('Douglas rachford needs exactly 2 functions')
    end
    
end


function [sol, s] = douglas_rachford_algorithm(Fp, sol, s, param)
%     s.x_n{1} = Fp{1}.prox_ad(2*sol-s.u_n,param.gamma);
%     s.u_n=s.u_n+param.lambda*(s.x_n{1}{1}-sol);
%     s.x_n{2} =Fp{2}.prox_ad(s.u_n,param.gamma);
%     sol = s.x_n{2}{1};
    s.x_n{1} = Fp{1}.prox(2*sol-s.u_n,param.gamma);
    s.u_n=s.u_n+param.lambda*(s.x_n{1}-sol);
    s.x_n{2} =Fp{2}.prox(s.u_n,param.gamma);
    sol = s.x_n{2};
       
end
