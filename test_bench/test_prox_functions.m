function errors = test_prox_functions()

errors = 0;
    
errors = errors + test_prox_adjoint();
errors = errors + test_prox_sumg();



end



function errors = test_prox_adjoint()

    errors = 0;
    tol = 1e-10;
    N = 10;
    
    x = 3*randn(N,1);
    
    param.verbose = 0;
    f.prox = @(x,T) prox_l1(x,T,param);
    f2.prox = @(x,T) prox_adjoint(x,T,f);
    T = rand;
    p1 = f.prox(x,T);
    p2 = prox_adjoint(x,T,f2);
    
    if norm(p1-p2)/norm(p1)<tol
        fprintf('  Test prox adjoint OK\n')
    else
        fprintf('  Test prox adjoint Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p1-p2)/norm(p1)
        errors= errors +1;
    end


end



function errors = test_prox_sumg()
    gsp_reset_seed();
    errors = 0;
    tol = 1e-4;
    N = 10;
    
    x0 = 3*randn(N,1);
    
    param.verbose = 0;
    f1.prox = @(x,T) prox_l1(x,T,param);
    f1.eval = @(x) norm(x,1);
    M = double( rand < 0.5 );
    f2.prox = @(x,T) (x - M.*x) + x0;
    f2.eval = @(x) eps;
    

    
    T = rand;

    f3.grad =  @(x) 1/T*(x-x0);
    f3.eval = @(x) 0.5/T*norm(x-x0);
    f3.beta = 1/T;
    
    param.G = {f1,f2};
    param.maxit = 100;
    param.tol = 0;
    param.verbose = 0;
    p1 = prox_sumg(x0,T,param);

    param.algo = 'generalized_forward_backward'; 

    p2 = solvep(x0,{f1,f2,f3},param);
    
    if norm(p1-p2)/norm(p1)<tol
        fprintf('  Test sumg OK\n')
    else
        fprintf('  Test sumg Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p1-p2)/norm(p1)
        errors= errors +1;
    end


end