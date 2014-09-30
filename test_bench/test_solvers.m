function [ errors ] = test_solvers( )
%TEST_SOLVERS test differents solvers
errors=0;

errors=errors+test_fwbw_simple();
errors=errors+test_dr_simple();
errors=errors+test_dr_complex();
errors=errors+test_admm_simple();
errors=errors+test_ppxa_simple();
errors=errors+test_ppxa_complex();



errors=errors+test_admm_complex();
errors=errors+test_admm_complex2();
errors=errors+test_admm_complex3();
errors=errors+test_admm_complex4();
errors=errors+test_admm_complex5();






end


function [errors]=test_fwbw_simple()
    errors = 0;
    
    N =100;
    
    y = 3*rand(N,1);
    x0 = rand(N,1);
    
    f1.eval = @(x) 1/2*norm(x-y)^2;
    paraml2.y = y;
    paraml2.verbose = 0;
    f1.prox = @(x,T) prox_l2(x,0.5*T,paraml2);
    f1.grad = @(x) x-y;
    
    paraml1.verbose = 0;
    f2.eval = @(x) norm(x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    
    param.tol = 100*eps;
    param.verbose = 0;
    
    p2 = forward_backward(x0, f2,f1, param);
    p3 = prox_l1(y,1,paraml1);
    
    if norm(p3-p2)/norm(p3)<1e-6
        fprintf('  Test fwbw simple 1 OK\n')
    else
        fprintf('  Test fwbw simple Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end

function [errors]=test_dr_simple()
    errors = 0;
    
    N =100;
    
    y = 3*rand(N,1);
    x0 = rand(N,1);
    
    f1.eval = @(x) 1/2*norm(x-y)^2;
    paraml2.y = y;
    paraml2.verbose = 0;
    f1.prox = @(x,T) prox_l2(x,0.5*T,paraml2);
    f1.grad = @(x) x-y;
    
    paraml1.verbose = 0;
    f2.eval = @(x) norm(x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    
    param.tol = 100*eps;
    param.verbose = 0;
    param.gamma = 1;
    
    
    p2 = douglas_rachford(x0,f1,f2, param);
    p3 = prox_l1(y,1,paraml1);
    
    if norm(p3-p2)/norm(p3)<1e-6
        fprintf('  Test dr simple 1 OK\n')
    else
        fprintf('  Test dr simple Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end

function [errors]=test_admm_simple()
    errors = 0;
    
    N =100;
    
    y = 3*rand(N,1);
    x0 = rand(N,1);
    
    f1.eval = @(x) 1/2*norm(x-y)^2;
    paraml2.y = y;
    paraml2.verbose = 0;
    f1.prox = @(x,T) prox_l2(x,0.5*T,paraml2);
    f1.grad = @(x) x-y;
    
    paraml1.verbose = 0;
    f2.eval = @(x) norm(x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    
    param.tol = 100*eps;
    param.verbose = 0;
    param.gamma = 1;
    
    
    p2 = admm(x0,f1,f2, param);
    p3 = prox_l1(y,1,paraml1);
    
    if norm(p3-p2)/norm(p3)<1e-6
        fprintf('  Test admm simple 1 OK\n')
    else
        fprintf('  Test admm simple Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end


function [errors]=test_ppxa_simple()
    errors = 0;
    
    N =100;
    
    y = 3*rand(N,1);
    x0 = rand(N,1);
    
    f1.eval = @(x) 1/2*norm(x-y)^2;
    paraml2.y = y;
    paraml2.verbose = 0;
    f1.prox = @(x,T) prox_l2(x,0.5*T,paraml2);
    f1.grad = @(x) x-y;
    
    paraml1.verbose = 0;
    f2.eval = @(x) norm(x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    
    param.tol = 100*eps;
    param.verbose = 0;
    param.gamma = 1;
    
    
    p2 = ppxa(x0,{f1,f2}, param);
    p3 = prox_l1(y,1,paraml1);
    
    if norm(p3-p2)/norm(p3)<1e-6
        fprintf('  Test ppxa simple 1 OK\n')
    else
        fprintf('  Test ppxa simple Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end


function [errors]=test_ppxa_complex()
    errors = 0;
    
    N =100;
    
    y = 3*rand(N,1);
    x0 = rand(N,1);
    
    
    f1.eval = @(x) 1/2*norm(x-y)^2;
    paraml2.y = y;
    paraml2.verbose = 0;
    f1.prox = @(x,T) prox_l2(x,0.5*T,paraml2);
    f1.grad = @(x) x-y;
    
    A = rand(N);
    
    paraml1.verbose = 0;
    paraml1.A = @(x) A*x;
    paraml1.At = @(x) A'*x;
    paraml1.nu = norm(A)^2;
    paraml1.tight = 0;
    f2.eval = @(x) norm(A*x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    
    param.tol = 100*eps;
    param.verbose = 0;
    param.gamma = 1;
    
    
    p2 = ppxa(x0,{f1,f2}, param);
    p3 = douglas_rachford(x0,f1,f2,param);
    
    if norm(p3-p2)/norm(p3)<1e-6
        fprintf('  Test ppxa complex 1 OK\n')
    else
        fprintf('  Test ppxa complex Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end


function [errors]=test_dr_complex()
    errors = 0;
    
    N =100;
    
    y = 3*rand(N,1);
    x0 = rand(N,1);
    
    
    f1.eval = @(x) 1/2*norm(x-y)^2;
    paraml2.y = y;
    paraml2.verbose = 0;
    f1.prox = @(x,T) prox_l2(x,0.5*T,paraml2);
    f1.grad = @(x) x-y;
    
    A = rand(N);
    
    paraml1.verbose = 0;
    paraml1.A = @(x) A*x;
    paraml1.At = @(x) A'*x;
    paraml1.nu = norm(A)^2;
    paraml1.tight = 0;
    f2.eval = @(x) norm(A*x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    
    param.tol = 100*eps;
    param.verbose = 0;
    param.gamma = 1;
    
    
    p2 = forward_backward(x0,f2,f1, param);
    p3 = douglas_rachford(x0,f1,f2,param);
    
    if norm(p3-p2)/norm(p3)<1e-6
        fprintf('  Test dr complex 1 OK\n')
    else
        fprintf('  Test dr complex Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end

function [errors]=test_admm_complex()
    errors = 0;
    
    N =100;
    
    y = 3*rand(N,1);
    x0 = rand(N,1);
    
    f1.eval = @(x) 1/2*norm(x-y)^2;
    f1.prox = @(x,T) (sqrt(N) * ifft(x)+y)/(1+T);
  
    paraml1.verbose = 0;
    f2.eval = @(x) norm(x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    
    param.tol = 100*eps;
    param.verbose = 0;
    param.gamma = 1;
    param.L =@(x) 1/sqrt(N) * fft(x);
    
    p2 = admm(x0,f1,f2, param);
    
    paraml1.A = @(x) 1/sqrt(N) * fft(x);
    paraml1.At = @(x) sqrt(N) * ifft(x);
    p3 = prox_l1(y,1,paraml1);
    
    if norm(p3-p2)/norm(p3)<1e-6
        fprintf('  Test admm complex 1 OK\n')
    else
        fprintf('  Test admm complex Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end



function [errors]=test_admm_complex2()
    errors = 0;
    
    N =100;
    
    y = 3*rand(N,1);
    x0 = rand(N,1);
    
    f1.eval = @(x) 1/2*norm(x-y)^2;
    
    paraml2.verbose = 0;
    paraml2.A = @(x) 1/sqrt(N) * fft(x);
    paraml2.At = @(x) sqrt(N) * ifft(x);
    paraml2.tight = 1;
    f1.prox = @(x,T) reverse_prox(x,0.5*T,paraml2,y);
  
    paraml1.verbose = 0;
    f2.eval = @(x) norm(x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    
    param.tol = 100*eps;
    param.verbose = 0;
    param.gamma = 1;
    param.L =@(x) 1/sqrt(N) * fft(x);
    
    p2 = admm(x0,f1,f2, param);
    
    paraml1.A = @(x) 1/sqrt(N) * fft(x);
    paraml1.At = @(x) sqrt(N) * ifft(x);
    p3 = prox_l1(y,1,paraml1);
    
    if norm(p3-p2)/norm(p3)<1e-6
        fprintf('  Test admm complex 2 OK\n')
    else
        fprintf('  Test admm complex 2 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end

function sol = reverse_prox(x,T,paraml2,z)
    paraml2.y = x;
    sol = prox_l2(z,T,paraml2);
end



function [errors]=test_admm_complex3()
    errors = 0;
    
    N =100;
    
    y = 3*rand(N,1);
    x0 = rand(N,1);
    
    f1.eval = @(x) 1/2*norm(x-y)^2;
    
    paraml2.verbose = 0;
    paraml2.A = @(x) 1/sqrt(N) * fft(x);
    paraml2.At = @(x) sqrt(N) * ifft(x);
    paraml2.nu = 10;
    paraml2.tight = 0;
    f1.prox = @(x,T) reverse_prox(x,0.5*T,paraml2,y);
  
    paraml1.verbose = 0;
    f2.eval = @(x) norm(x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    
    param.tol = 100*eps;
    param.verbose = 0;
    param.gamma = 1;
    param.L =@(x) 1/sqrt(N) * fft(x);
    
    p2 = admm(x0,f1,f2, param);
    
    paraml1.A = @(x) 1/sqrt(N) * fft(x);
    paraml1.At = @(x) sqrt(N) * ifft(x);
    p3 = prox_l1(y,1,paraml1);
    
    if norm(p3-p2)/norm(p3)<1e-6
        fprintf('  Test admm complex 3 OK\n')
    else
        fprintf('  Test admm complex 3 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end



function [errors]=test_admm_complex4()
    errors = 0;
    
    N =100;
    
    y = 3*rand(N,1);
    x0 = rand(N,1);
    
    a=25;
    M=50;
    g=gabwin({'gauss',1},a,M,N);
    gt=gabtight(g,a,M,N);
    gt=gt/norm(gt);
    Fr=frame('dgt',gt,a,M);
    F=frsynmatrix(Fr,N);
    [~,B]=gabframebounds(gt,a,M,N);



    
    
    A = @(x) F'*x;%  @(x) gradient_op1d(x);
    At = @(x) F*x;% @(x) -div_op1d(x);
    
    f1.eval = @(x) 1/2*norm(x-y)^2;
    paraml2.verbose = 0;
    paraml2.A = A;
    paraml2.At = At;
    paraml2.nu = B;
    paraml2.tight = 0;
    f1.prox = @(x,T) reverse_prox(x,0.5*T,paraml2,y);
  
    paraml1.verbose = 0;
    f2.eval = @(x) norm(x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    param.maxit = 1000;
    param.tol = 0;
    param.verbose = 0;
    param.L =A;
    
    p2 = admm(x0,f1,f2, param);
    
    paraml1.verbose = 0;
    paraml1.maxit = 1000;
    paraml1.tight = 0;
    paraml1.tol = 0;
    paraml1.A = A;
    paraml1.At = At;
    paraml1.nu = B;
    p3 = prox_l1(y,1,paraml1);
    
    if norm(p3-p2)/norm(p3)<2e-3
        fprintf('  Test admm complex 4 OK\n')
    else
        fprintf('  Test admm complex 4 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end




function [errors]=test_admm_complex5()
    errors = 0;
    
    N =100;
    
    y = 3*rand(N,1);
    x0 = rand(N,1);

    
    
    A =   @(x) gradient_op1d(x);
    At = @(x) -div_op1d(x);
    
    f1.eval = @(x) 1/2*norm(x-y)^2;
    paraml2.verbose = 0;
    paraml2.A = A;
    paraml2.At = At;
    paraml2.nu = 4;
    paraml2.tight = 0;
    f1.prox = @(x,T) reverse_prox(x,0.5*T,paraml2,y);
  
    paraml1.verbose = 0;
    f2.eval = @(x) norm(x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    param.maxit = 1000;
    param.tol = 0;
    param.verbose = 0;
    param.L =A;
    
    p2 = admm(x0,f1,f2, param);
    
    paraml1.verbose = 0;
    paraml1.maxit = 1000;
    paraml1.tol = 0;
    p3 = prox_tv1d(y,1,paraml1);
    
    if norm(p3-p2)/norm(p3)<2e-3
        fprintf('  Test admm complex 5 OK\n')
    else
        fprintf('  Test admm complex 5 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end
