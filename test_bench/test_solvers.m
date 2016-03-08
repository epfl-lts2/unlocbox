function [ errors ] = test_solvers( )
%TEST_SOLVERS test differents solvers
errors=0;
gsp_reset_seed(1);

errors=errors+test_chambolle_pock_simple();
errors=errors+test_chambolle_pock_complex();


errors = errors+ test_primal_dual();


errors=errors+test_all_solver();



errors=errors+test_pocs();
errors=errors+test_gradient_descent();

errors = errors+test_gefwbw_simple();
errors = errors+test_gefwbw_with_old();




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


errors = errors + test_sdmm_simple();
errors = errors + test_sdmm_complex();
errors = errors + test_sdmm_complex2();


errors = errors+ test_fb_based_primal_dual();
errors = errors+ test_fb_based_primal_dual2();
errors = errors+ test_fb_based_primal_dual3();

errors = errors+ test_fb_based_primal_dual4();

errors = errors+ test_fb_based_primal_dual_fista();


errors=errors + test_verbosity();


end

function errors = test_pocs()
    errors = 0;
    x0 = [0 0]';

    x1 = [3 3]';
    x2 = [-3 3]';

    paramb21.verbose = 0;
    paramb21.y = x1;
    paramb21.epsilon = 4;
    f1.prox = @(x,T) proj_b2(x,T,paramb21);
    f1.eval = @(x) max(norm(x-x1)-4,0);

    paramb22.verbose = 0;
    paramb22.y = x2;
    paramb22.epsilon = 4;
    f2.prox =@(x,T) proj_b2(x,T,paramb22);
    f2.eval = @(x) max(norm(x-x2)-4,0);
    
    param.tol = 100 *eps;
    param.maxit = 1000;
    param.verbose = 0;
    sol = pocs(x0,{f1,f2},param);
    solr = [0, 3-sqrt(7)]';

    if norm(sol-solr)/norm(solr)<1e-9
        fprintf('  Test pocs OK\n')
    else
        fprintf('  Test pocs Pas OK!!!!!!!!!!!!!!!!\n')
        norm(sol-solr)/norm(solr)
        errors= errors +1;
    end


end

function errors = test_gradient_descent()
    errors = 0;
    N = 20;
    M = 10;
    A = rand(N,M);
    y = rand(N,1);
    x0 = zeros(M,1);
    sol = pinv(A) * y;
    
    f.eval = @(x) norm(A*x-y);
    f.grad = @(x) 2*A'*(A*x-y);
    f.beta = 2*norm(A)^2;
    
    param.tol = 10 * eps;
    param.maxit = 2001;
    param.verbose = 0;
    sol2 = gradient_descent(x0,f,param);

    if norm(sol2-sol)/norm(sol)<2e-3
        fprintf('  Test gradient descent OK\n')
    else
        fprintf('  Test gradient descent Pas OK!!!!!!!!!!!!!!!!\n')
        norm(sol2-sol)/norm(sol)
        errors= errors +1;
    end
end


function [errors]=test_all_solver()
    errors = 0;
    N = 50;
    img = cameraman();
    img = img(100:100+N,100:100+N);
    
    M = rand(size(img));
    M = M>0.5;
    
    A = @(x) M.*x;
    At = A;
    
    y = A(img) + 0.03 * randn(size(img));
    
    lambda = 0.03;

    x0 = y;
    H =@(x) 1/N* dct2(x);
    Ht = @(x) N * idct2(x);
    
    f1.eval = @(x) 1/2*norm(A(x)-y)^2;
    f1.grad = @(x) At(A(x)-y);
    paraml2.A = A;
    paraml2.At = At;
    paraml2.y = y;
    paraml2.verbose = 0;
    paraml2.tight = 0;
    f1.prox = @(x,T) prox_l2(x,0.5*T,paraml2);
    f1.beta = 1;
    
    paraml1.verbose = 0;
    paraml1.A = H;
    paraml1.At = Ht;
    f2.eval = @(x) lambda*norm(H(x),1);
    f2.prox = @(x,T) prox_l1(x,T*lambda,paraml1);
    
    param.tol = 100*eps;
    param.verbose = 0;
    param.maxit = 1000;
    
    f3.prox = @(x,T) x;
    f3.eval = @(x) eps;
    
    p1 = forward_backward(x0,f2,f1,param);

    p2 = generalized_forward_backward(x0, {f2},f1, param);
    param.gamma = 1;

    p3 = douglas_rachford(x0,f2,f1,param);
    p5 = admm(x0,f2,f1, param);

    f1.x0 = x0;

    param.maxit = 2000;
    p4 = ppxa(x0, {f2,f1}, param);
    
    paraml12.verbose = 0;
    f22.eval = @(x) lambda*norm(H(x),1);
    f22.prox = @(x,T) prox_l1(x,T*lambda,paraml12);
    f22.x0 = H(x0);
    f22.L = H;
    f22.Lt = Ht;
    f22.norm_L = 1;
    
    p6 = sdmm({f22,f1}, param);
    %p7 = chambolle_pock(x0, f2,f1, param);

    p8 = fb_based_primal_dual(x0,f1, f22, f3, param);
    
    
    if norm(p2-p1)/norm(p1)<1e-2
        fprintf('  Test all gefwbw OK\n')
    else
        fprintf('  Test all gefwbw Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p2-p1)/norm(p1)
        errors= errors +1;
    end
    
    if norm(p3-p1)/norm(p1)<1e-2
        fprintf('  Test all dg OK\n')
    else
        fprintf('  Test all dg Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p1)/norm(p1)
        errors= errors +1;
    end
    
    if norm(p4-p1)/norm(p1)<1e-2
        fprintf('  Test all ppxa OK\n')
    else
        fprintf('  Test all ppxa Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p4-p1)/norm(p1)
        errors= errors +1;
    end

    if norm(p5-p1)/norm(p1)<1e-2
        fprintf('  Test all admm 1 OK\n')
    else
        fprintf('  Test all admm 1 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p5-p1)/norm(p1)
        errors= errors +1;
    end
    
    
    if norm(p6-p1)/norm(p1)<1e-2
        fprintf('  Test all sdmm 1 OK\n')
    else
        fprintf('  Test all sdmm 1 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p6-p1)/norm(p1)
        errors= errors +1;
    end
    
%     if norm(p7-p1)/norm(p1)<1e-5
%         fprintf('  Test all chambolle pock 1 OK\n')
%     else
%         fprintf('  Test all chambolle pock 1 Pas OK!!!!!!!!!!!!!!!!\n')
%         norm(p7-p1)/norm(p1)
%         errors= errors +1;
%     end

    if norm(p8-p1)/norm(p1)<1e-2
        fprintf('  Test all fb_based_primal_dual 1 OK\n')
    else
        fprintf('  Test all fb_based_primal_dual 1 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p8-p1)/norm(p1)
        errors= errors +1;
    end
end

function [errors]=test_gefwbw_simple()
    errors = 0;
    
    N =100;
    
    y = 3*rand(N,1);
    x0 = rand(N,1);
    
    f1.eval = @(x) 1/2*norm(x-y)^2;
    f1.grad = @(x) x-y;
    f1.beta = 1;
    
    paraml1.verbose = 0;
    f2.eval = @(x) norm(x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    
    param.tol = 100*eps;
    param.verbose = 0;
    
    p2 = generalized_forward_backward(x0, {f2},f1, param);
    p3 = prox_l1(y,1,paraml1);
    
    if norm(p3-p2)/norm(p3)<1e-6
        fprintf('  Test gefwbw simple 1 OK\n')
    else
        fprintf('  Test gefwbw simple 1 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end

function [errors]=test_gefwbw_with_old()
    errors = 0;
    
    N =10;
    
    y = 3*rand(N,1);
    x0 = rand(N,1);
    A = rand(N);
    
    f1.eval = @(x) norm(A*x-y)^2;
    f1.grad = @(x) 2*A'*(A*x-y);
    f1.beta = 2*norm(A)^2;
    
    paraml1.verbose = 0;
    f2.eval = @(x) norm(x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    
    f3 = f2;
    f4 = f2;
    param.maxit = 3000;
    param.tol = 100*eps;
    param.verbose = 0;
    p3 = solvep(x0, {f2,f3,f1,f4}, param);




    f5.eval = @(x) norm(3*x,1);
    f5.prox = @(x,T) prox_l1(x,3*T,paraml1);
    p4 = solvep(x0, {f5,f1}, param);
    
    f6.eval = @(x) norm(2*x,1);
    f6.prox = @(x,T) prox_l1(x,2*T,paraml1);
    param.algo = 'FB_BASED_PRIMAL_DUAL';
    p5 = solvep(x0, {f1,f4,f6}, param);
    
    param.gamma = 1/(2*norm(A)^2);
    
    p2 = generalized_forward_backward_old(x0, {f2,f3,f4},f1, param);
    
    if norm(p3-p2)/norm(p3)<1e-6
        fprintf('  Test gefwbw with old 1 OK\n')
    else
        fprintf('  Test gefwbw with old 1 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
    if norm(p3-p4)/norm(p3)<1e-5
        fprintf('  Test gefwbw with FISTA 1 OK\n')
    else
        fprintf('  Test gefwbw with FISTA 1 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p4)/norm(p3)
        errors= errors +1;
    end
    
    if norm(p5-p4)/norm(p4)<1e-5
        fprintf('  Test primal dual with FISTA 1 OK\n')
    else
        fprintf('  Test primal dual with FISTA 1 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p5-p4)/norm(p4)
        errors= errors +1;
    end
    
end

function [errors]=test_fwbw_simple()
    errors = 0;
    
    N =100;
    
    y = 3*rand(N,1);
    x0 = rand(N,1);
    
    f1.eval = @(x) 1/2*norm(x-y)^2;
    f1.grad = @(x) x-y;
    f1.beta = 1;
    
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
        fprintf('  Test fwbw simple 1 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    N = 10;
    
    A = dct(eye(N));
    y = N*rand(N,1);
    x0 = rand(N,1);

    f1.eval = @(x) 1/2*norm(A*x-y)^2;
    f1.grad = @(x) A'*(A*x-y);
    f1.beta = norm(A)^2;
    param.method = 'FISTA';
    param.maxit = 1000;
    param.tol = eps;

    p5 = forward_backward(x0, f2,f1,param);

    param.method = 'ISTA';
    p4 = forward_backward(x0, f2,f1,param);

    
    if norm(p4-p5)/norm(p4)<1e-5
        fprintf('  Test fwbw 2 OK\n')
    else
        fprintf('  Test fwbw 2 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p4-p5)/norm(p4)
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
    f1.beta = 1;
    
    paraml1.verbose = 0;
    f2.eval = @(x) norm(x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    
    param.tol = 100*eps;
    param.verbose = 0;
    param.gamma = 1;
    
    
    [p2,infos] = admm(x0,f1,f2, param);
    p3 = prox_l1(y,1,paraml1);
    
    if norm(p3-p2)/norm(p3)<1e-6
        fprintf('  Test admm simple 1 OK\n')
    else
        fprintf('  Test admm simple Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    if strcmp(infos.crit,'REL_NORM_PRIMAL_DUAL')
        fprintf('  Test admm simple 2 OK\n')
    else
        fprintf('  Test admm simple 2 Pas OK!!!!!!!!!!!!!!!!\n')
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
    param.maxit = 1000;
    
    
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
    f1.beta = 1;
    
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
    
    L = @(x) 1/sqrt(N) * fft(x);
    
    f1.eval = @(x) 1/2*norm(x-y)^2;
    f1.proxL = @(x,T) (sqrt(N) * ifft(x)+y)/(1+T);
    
    paraml1.verbose = 0;
    
    f2.eval = @(x) norm(x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    f2.L = L;
   
    param.tol = 100*eps;
    param.verbose = 0;
    param.gamma = 1;

    
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
    f1.proxL = @(x,T) reverse_prox(x,0.5*T,paraml2,y);
    f1.prox = @(x,T) prox_l2(x,T);
    
    paraml1.verbose = 0;
    f2.eval = @(x) norm(x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    f2.L = @(x) 1/sqrt(N) * fft(x);
    
    param.tol = 100*eps;
    param.verbose = 0;
    param.gamma = 1;
    
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
    f1.proxL = @(x,T) reverse_prox(x,0.5*T,paraml2,y);
  
    paraml1.verbose = 0;
    f2.eval = @(x) norm(x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    f2.L = @(x) 1/sqrt(N) * fft(x);
    
   
    param.tol = 100*eps;
    param.verbose = 0;
    param.gamma = 1;
    
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
    f1.proxL = @(x,T) reverse_prox(x,0.5*T,paraml2,y);
  
    paraml1.verbose = 0;
    f2.eval = @(x) norm(x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    f2.L = A;
    param.maxit = 1000;
    param.tol = 0;
    param.verbose = 0;
    
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
    f1.proxL = @(x,T) reverse_prox(x,0.5*T,paraml2,y);
    
    paraml1.verbose = 0;
    f2.eval = @(x) norm(x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    f2.L =A;

    param.maxit = 1000;
    param.tol = 0;
    param.verbose = 0;
    
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


function [errors]=test_sdmm_simple()
    errors = 0;
    
    N =100;
    
    y = 3*rand(N,1);
    x0 = rand(N,1);
    
    f1.eval = @(x) 1/2*norm(x-y)^2;
    paraml2.y = y;
    paraml2.verbose = 0;
    f1.prox = @(x,T) prox_l2(x,0.5*T,paraml2);
    f1.grad = @(x) x-y;
    f1.beta = 1;
    f1.x0 = x0;
    
    paraml1.verbose = 0;
    f2.eval = @(x) norm(x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    f2.x0 = x0;
    
    param.tol = 1e-10;
    param.verbose = 0;
    param.gamma = 1;
    
    
    [p2,infos] = sdmm({f1,f2}, param);
    p3 = prox_l1(y,1,paraml1);
    
    if norm(p3-p2)/norm(p3)<1e-6
        fprintf('  Test sdmm simple 1 OK\n')
    else
        fprintf('  Test sdmm simple Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
end

function [errors]=test_sdmm_complex()
    errors = 0;
    
    N =100;
    
    y = 3*rand(N,1);
    x0 = rand(N,1);
    
    paraml2.verbose = 0;
    paraml2.y = y;
    f1.eval = @(x) 1/2*norm(x-y)^2;
    f1.prox = @(x,T) prox_l2(x,0.5*T,paraml2);
    f1.x0 = x0;
  
    paraml1.verbose = 0;
    f2.eval = @(x) norm(x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    f2.L = @(x) 1/sqrt(N) * fft(x);
    f2.Lt = @(x) sqrt(N) * ifft(x);
    f2.x0 = zeros(size(x0));
    
    param.tol = 100*eps;
    param.verbose = 0;
    param.gamma = 1;
    param.Qinv = @(x) 1/2 *x;
    
    p2 = sdmm({f1,f2}, param);
    
    paraml1.A = @(x) 1/sqrt(N) * fft(x);
    paraml1.At = @(x) sqrt(N) * ifft(x);
    p3 = prox_l1(y,1,paraml1);
    
    if norm(p3-p2)/norm(p3)<1e-6
        fprintf('  Test sdmm complex 1 OK\n')
    else
        fprintf('  Test sdmm complex Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end




function [errors]=test_sdmm_complex2()
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



    
    
    A = @(x) F'*x;
    At = @(x) F*x;
    
    f1.eval = @(x) 1/2*norm(x-y)^2;
    paraml2.verbose = 0;
    paraml2.y = y;
    f1.prox = @(x,T) prox_l2(x,0.5*T,paraml2);
    f1.x0 = zeros(size(y));
  
    paraml1.verbose = 0;
    f2.eval = @(x) norm(x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    f2.L = F';
    f2.Lt = F;
    f2.x0 = zeros(M/a*N,1);
    param.maxit = 1000;
    param.tol = 0;
    param.verbose = 0;
    param.Qinv = @(x) 1/(1+B) * x;
    
    p2 = sdmm({f1,f2}, param);
    
    paraml1.verbose = 0;
    paraml1.maxit = 1000;
    paraml1.tight = 0;
    paraml1.tol = 0;
    paraml1.A = A;
    paraml1.At = At;
    paraml1.nu = B;
    p3 = prox_l1(y,1,paraml1);
    
    if norm(p3-p2)/norm(p3)<2e-3
        fprintf('  Test sdmm complex 2 OK\n')
    else
        fprintf('  Test sdmm complex 2 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end

function [errors]=test_chambolle_pock_complex()
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
    paraml2.y = y;
    f1.prox = @(x,T) prox_l2(x,0.5*T,paraml2);

  
    paraml1.verbose = 0;
    f2.eval = @(x) norm(x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    f2.L = A;
    f2.Lt = At;
    f2.norm_L = B;
    param.maxit = 1000;
    param.tol = 1000 *eps;
    param.verbose = 0;

    
    p2 = chambolle_pock(x0,f1,f2, param);
    
    paraml1.verbose = 0;
    paraml1.maxit = 1000;
    paraml1.tight = 0;
    paraml1.tol = 0;
    paraml1.A = A;
    paraml1.At = At;
    paraml1.nu = B;
    p3 = prox_l1(y,1,paraml1);
    
    if norm(p3-p2)/norm(p3)<2e-3
        fprintf('  Test chambolle pock complex OK\n')
    else
        fprintf('  Test chambolle pock complex Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end



function [errors] = test_chambolle_pock_simple()
    errors = 0;
    
    N =5;
    
    y = 3*rand(N,1);
    x0 = rand(N,1);
    
    f1.eval = @(x) 1/2*norm(x-y)^2;
    paraml2.y = y;
    paraml2.verbose = 0;
    f1.prox = @(x,T) prox_l2(x,0.5*T,paraml2);
    f1.grad = @(x) x-y;
    f1.beta = 1;
    
    paraml1.verbose = 0;
    f2.eval = @(x) norm(x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    
    param.tol = 100*eps;
    param.verbose = 0;
    param.gamma = 1;
    
    
    [p2,infos] = chambolle_pock(x0,f1,f2, param);
    p3 = solvep(x0, {f1,f2},param);
    
    if norm(p3-p2)/norm(p3)<1e-6
        fprintf('  Test chambolle pock  simple 1 OK\n')
    else
        fprintf('  Test chambolle pock  simple Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
    if strcmp(infos.crit,'REL_NORM_PRIMAL_DUAL')
        fprintf('  Test chambolle pock simple 2 OK\n')
    else
        fprintf('  Test chambolle pock simple 2 Pas OK!!!!!!!!!!!!!!!!\n')
        errors= errors +1;
    end
end



function [errors] = test_fb_based_primal_dual()
    errors = 0;
    
    N =10;
    
    y = 3*rand(N,1);
    x0 = rand(N,1);
    
    f1.eval = @(x) 1/2*norm(x-y)^2;
    paraml2.y = y;
    paraml2.verbose = 0;
    f1.prox = @(x,T) prox_l2(x,0.5*T,paraml2);
    f1.grad = @(x) x-y;
    f1.beta = 1;
    
    paraml1.verbose = 0;
    f2.eval = @(x) norm(x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    
    param.tol = 100*eps;
    param.verbose = 0;
    param.gamma = 1;

%     
    f3.prox = @(x,T) x;
    f3.eval = @(x) eps;
    
    [p2] = fb_based_primal_dual(x0,f1,f2,f3, param);
    p3 = prox_l1(y,1,paraml1);

    
    if norm(p3-p2)/norm(p3)<1e-6
        fprintf('  Test fb_based_primal_dual simple 1 OK\n')
    else
        fprintf('  Test fb_based_primal_dual  simple Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end




function [errors] = test_fb_based_primal_dual2()
    errors = 0;
    
    N =10;
    
    y = 3*rand(N,1);
    x0 = rand(N,1);
    
    f1.eval = @(x) 1/2*norm(x-y)^2;
    paraml2.y = y;
    paraml2.verbose = 0;
    f1.prox = @(x,T) prox_l2(x,0.5*T,paraml2);
    f1.grad = @(x) x-y;
    f1.beta = 1;
    
    paraml1.verbose = 0;
    f2.eval = @(x) norm(x,1);
    f2.prox = @(x,T) prox_l1(x,T,paraml1);
    
    param.tol = 100*eps;
    param.verbose = 0;
    param.gamma = 0.1;
    param.maxit = 1000;
     
    paramb2.verbose = 0;
    paramb2.y = 5*rand(N,1);
    paramb2.epsilon = 2;
    f3.prox = @(x,T) proj_b2(x,T,paramb2);
    f3.eval = @(x) eps;

    
    p2 = fb_based_primal_dual(x0,f1,f2,f3, param);

    p3 = ppxa(x0,{f1,f2,f3},param);
    
    if norm(p3-p2)/norm(p3)<1e-6
        fprintf('  Test fb_based_primal_dual 2 OK\n')
    else
        fprintf('  Test fb_based_primal_dual  2 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end


function [errors] = test_fb_based_primal_dual3()
    errors = 0;
    
    N =10;
    
    y = 3*rand(N,1);
    x0 = rand(N,1);
    A =@(x) 1/sqrt(N) * dct(x);
    At = @(x) sqrt(N) * idct(x);
    f1.eval = @(x) 1/2*norm(x-y)^2;
    paraml2.y = y;
    paraml2.verbose = 0;
    f1.prox = @(x,T) prox_l2(x,0.5*T,paraml2);
    f1.grad = @(x) x-y;
    f1.beta = 1;
    
    paraml11.verbose = 0;
    f21.eval = @(x) norm(x,1);
    f21.prox = @(x,T) prox_l1(x,T,paraml11);    
    f21.L = A;
    f21.Lt = At;
    f22.norm_L = 1;
    
    paraml12.verbose = 0;
    paraml12.A = A;
    paraml12.At = At;
    f22.eval = @(x) norm(A(x),1);
    f22.prox = @(x,T) prox_l1(x,T,paraml12);
    
    param.tol = 100*eps;
    param.verbose = 0;
    param.gamma = 0.1;
    param.maxit = 1000;
     

    f3.prox = @(x,T) x;
    f3.eval = @(x) eps;
    
    p2 = fb_based_primal_dual(x0,f1,f21,f3, param);

    p3 = f22.prox(y,1);
    
    if norm(p3-p2)/norm(p3)<1e-6
        fprintf('  Test fb_based_primal_dual 3 OK\n')
    else
        fprintf('  Test fb_based_primal_dual  3 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end



function [errors] = test_fb_based_primal_dual4()
    errors = 0;
    
    N =10;
    
    y = 3*rand(N,1);
    x0 = rand(N,1);
    A =@(x) 1/sqrt(N) * dct(x);
    At = @(x) sqrt(N) * idct(x);
    f1.eval = @(x) 1/2*norm(x-y)^2;
    paraml2.y = y;
    paraml2.verbose = 0;
    f1.prox = @(x,T) prox_l2(x,0.5*T,paraml2);
    f1.grad = @(x) x-y;
    f1.beta = 1;
    
    paraml11.verbose = 0;
    f21.eval = @(x) norm(x,1);
    f21.prox = @(x,T) prox_l1(x,T,paraml11);    
    f21.L = A;
    f21.Lt = At;
    f21.norm_L = 1;
    
    paraml12.verbose = 0;
    paraml12.A = A;
    paraml12.At = At;
    f22.eval = @(x) norm(A(x),1);
    f22.prox = @(x,T) prox_l1(x,T,paraml12);
    
    param.tol = 100*eps;
    param.verbose = 0;
    param.gamma = 0.1;
    param.maxit = 1000;
     
    paramb2.verbose = 0;
    paramb2.y = 5*rand(N,1);
    paramb2.epsilon = 2;
    f3.prox = @(x,T) proj_b2(x,T,paramb2);
    f3.eval = @(x) eps;
    
    p2 = fb_based_primal_dual(x0,f1,f21,f3, param);

	p3 = ppxa(x0,{f1,f22,f3},param);
    
    if norm(p3-p2)/norm(p3)<1e-6
        fprintf('  Test fb_based_primal_dual 4 OK\n')
    else
        fprintf('  Test fb_based_primal_dual  4 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end



function [errors] = test_fb_based_primal_dual_fista()
    errors = 0;
    
    N =10;
    
    y = 3*rand(N,1);
    x0 = rand(N,1);
    A =@(x) 1/sqrt(N) * dct(x);
    At = @(x) sqrt(N) * idct(x);
    f1.eval = @(x) 1/2*norm(x-y)^2;
    paraml2.y = y;
    paraml2.verbose = 0;
    f1.prox = @(x,T) prox_l2(x,0.5*T,paraml2);
    f1.grad = @(x) x-y;
    f1.beta = 1;
    
    paraml11.verbose = 0;
    f21.eval = @(x) norm(x,1);
    f21.prox = @(x,T) prox_l1(x,T,paraml11);    
    f21.L = A;
    f21.Lt = At;
    f21.norm_L = 1;
    
    paraml12.verbose = 0;
    paraml12.A = A;
    paraml12.At = At;
    f22.eval = @(x) norm(A(x),1);
    f22.prox = @(x,T) prox_l1(x,T,paraml12);
    
    param.tol = 100*eps;
    param.verbose = 1;
    param.gamma = 0.1;
    param.maxit = 1000;
     
    paramb2.verbose = 0;
    paramb2.y = 5*rand(N,1);
    paramb2.epsilon = 2;
    f3.prox = @(x,T) proj_b2(x,T,paramb2);
    f3.eval = @(x) eps;
    param.method = 'FISTA';
    p2 = fb_based_primal_dual(x0,f1,f21,f3, param);
    param.method = 'ISTA';
	p3 = fb_based_primal_dual(x0,f1,f21,f3, param);
    
    if norm(p3-p2)/norm(p3)<1e-6
        fprintf('  Test fb_based_primal_dual fista OK\n')
    else
        fprintf('  Test fb_based_primal_dual  fista Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end


function errors = test_primal_dual()

    errors = 0;
    
    N =10;
    M = 6;
    y = 3*rand(N,1);
    x0 = rand(N,1);
    Mat = randn(M,N);
    A = @(x) Mat*x;
    At = @(x) Mat'*x;
    f1.eval = @(x) 1/2*norm(x-y)^2;
    paraml2.y = y;
    paraml2.verbose = 0;
    f1.prox = @(x,T) prox_l2(x,0.5*T,paraml2);
    f1.grad = @(x) x-y;
    f1.beta = 1;
    
    paraml11.verbose = 0;
    f21.eval = @(x) norm(x,1);
    f21.prox = @(x,T) prox_l1(x,T,paraml11);    
    f21.L = A;
    f21.Lt = At;
    f21.norm_L = norm(Mat)^2;
    
    paraml12.verbose = 0;
    f22.eval = @(x) norm(x,1);
    f22.prox = @(x,T) prox_l1(x,T,paraml12);
    
    param.tol = 100*eps;
    param.verbose = 1;
    param.maxit = 1000;
     
    paramb2.verbose = 0;
    paramb2.y = 5*rand(N,1);
    paramb2.epsilon = 2;
    f3.prox = @(x,T) proj_b2(x,T,paramb2);
    f3.eval = @(x) eps;

	p3 = fb_based_primal_dual(x0,f1,f21,f3, param);
    
    p2 = fbf_primal_dual(x0,f1,f21,f3, param);

    if norm(p3-p2)/norm(p3)<1e-6
        fprintf('  Test primal dual OK\n')
    else
        fprintf('  Test primal dual Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end


end


function errors = test_verbosity()
errors = 0;


x = randn(10,1);
y = randn(10,1);
A = rand(10);

paraml1.verbose = 0;
f1.prox = @(x,T) prox_l1(x,T,paraml1);
f1.eval = @(x) norm(x,1);

f2.prox = @(x,T) prox_l1(x,T,paraml1);
f2.eval = @(x) norm(x,1);

f3.grad = @(x,T) 2*A'*(A*x-y);
f3.eval = @(x) norm(A*x-y,2)^2;
f3.beta = 2*norm(A)^2;

param.maxit = 10;

try
    param.test_type = 'rel_norm_obj';
    param.verbose = 1;
    fb_based_primal_dual(x,f1,f2,f3,param);    
    param.verbose = 2;
    fb_based_primal_dual(x,f1,f2,f3,param);    
catch
    warning('Error in verbosity rel_norm_obj')
    errors = errors + 1;
end

try
    param.test_type = 'rel_norm_primal';
    param.verbose = 1;
    fb_based_primal_dual(x,f1,f2,f3,param);    
    param.verbose = 2;
    fb_based_primal_dual(x,f1,f2,f3,param);    
catch
    warning('Error in verbosity rel_norm_primal')
    errors = errors + 1;
end

try
    param.test_type = 'rel_norm_dual';
    param.verbose = 1;
    fb_based_primal_dual(x,f1,f2,f3,param);    
    param.verbose = 2;
    fb_based_primal_dual(x,f1,f2,f3,param);    
catch
    warning('Error in verbosity rel_norm_dual')
    errors = errors + 1;
end

try
    param.test_type = 'obj_increase';
    param.verbose = 1;
    fb_based_primal_dual(x,f1,f2,f3,param);    
    param.verbose = 2;
    fb_based_primal_dual(x,f1,f2,f3,param);    
catch
    warning('Error in verbosity obj_increase')
    errors = errors + 1;
end

try
    param.test_type = 'obj_threshold';
    param.verbose = 1;
    fb_based_primal_dual(x,f1,f2,f3,param);    
    param.verbose = 2;
    fb_based_primal_dual(x,f1,f2,f3,param);    
catch
    warning('Error in verbosity obj_threshold')
    errors = errors + 1;
end

param.debug_mode = 1;

try
    param.test_type = 'rel_norm_obj';
    param.verbose = 1;
    fb_based_primal_dual(x,f1,f2,f3,param);    
    param.verbose = 2;
    fb_based_primal_dual(x,f1,f2,f3,param);    
catch
    warning('Error in verbosity rel_norm_obj')
    errors = errors + 1;
end

try
    param.test_type = 'rel_norm_primal';
    param.verbose = 1;
    fb_based_primal_dual(x,f1,f2,f3,param);    
    param.verbose = 2;
    fb_based_primal_dual(x,f1,f2,f3,param);    
catch
    warning('Error in verbosity rel_norm_primal')
    errors = errors + 1;
end

try
    param.test_type = 'rel_norm_dual';
    param.verbose = 1;
    fb_based_primal_dual(x,f1,f2,f3,param);    
    param.verbose = 2;
    fb_based_primal_dual(x,f1,f2,f3,param);    
catch
    warning('Error in verbosity rel_norm_dual')
    errors = errors + 1;
end

try
    param.test_type = 'obj_increase';
    param.verbose = 1;
    fb_based_primal_dual(x,f1,f2,f3,param);    
    param.verbose = 2;
    fb_based_primal_dual(x,f1,f2,f3,param);    
catch
    warning('Error in verbosity obj_increase')
    errors = errors + 1;
end

try
    param.test_type = 'obj_threshold';
    param.verbose = 1;
    fb_based_primal_dual(x,f1,f2,f3,param);    
    param.verbose = 2;
    fb_based_primal_dual(x,f1,f2,f3,param);    
catch
    warning('Error in verbosity obj_threshold')
    errors = errors + 1;
end

end
