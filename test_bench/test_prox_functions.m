function errors = test_prox_functions()
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