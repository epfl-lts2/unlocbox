function [ errors ] = test_prox_l1( )
%TEST_PROX_L1 test prox_l1 and related function
errors=0;

errors=errors+test_simple(eps(100));
errors=errors+test_tight(eps(1000));

errors=errors+test_l1tv(eps(1000));
errors=errors+test_fax(eps(1000));

errors=errors+test_l1_y();




end

function [errors]=test_simple(tol)
    errors = 0;
    
    s = 2*rand(100);
    y = s;
    y(abs(y)<=1) = 0;
    y(abs(y)>0) = sign(y(abs(y)>0)) .* (abs(y(abs(y)>0)) - 1);
    
    param.verbose = 0;
    p2 = prox_l1(s,1,param);
    p3 = y; 
    
    if norm(p3-p2)/norm(p3)<tol
        fprintf('  Test prox_l1 simple 1 OK\n')
    else
        fprintf('  Test prox_l1 simple Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end



function [errors]=test_tight(tol)
    errors = 0;
    
    N = 100;
    
    s = 2*rand(N);
    y = 1/N*fft2(s);
    y(abs(y)<=1) = 0;
    y(abs(y)>0) = exp(1i*phase(y(abs(y)>0))) .* (abs(y(abs(y)>0)) - 1);
    
    param.verbose = 0;
    param.A = @(x) 1/N*fft2(x); 
    param.At = @(x) N*ifft2(x);
    
    p2 = prox_l1(s,1,param);
    p3 = N*ifft2(y); 
    
    if norm(p3-p2)/norm(p3)<tol
        fprintf('  Test prox_l1 tight 1 OK\n')
    else
        fprintf('  Test prox_l1 tight 1 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
    param.tight = 0;   
    p2 = prox_l1(s,1,param);

    
    if norm(p3-p2)/norm(p3)<tol
        fprintf('  Test prox_l1 tight 2 OK\n')
    else
        fprintf('  Test prox_l1 tight 2 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
    param.tight = 0;   
    param.nu = 5;
    param.tol = tol;
    p2 = prox_l1(s,1,param);

    
    if norm(p3-p2)/norm(p3)<100*tol
        fprintf('  Test prox_l1 tight 3 OK\n')
    else
        fprintf('  Test prox_l1 tight 3 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end


function [errors]=test_l1tv(tol)

    errors = 0;
    
    s = 10*rand(50,1);

    
    param.verbose = 0;
    param.A = @(x) gradient_op1d(x);
    param.At = @(x) -div_op1d(x);
    param.tight = 0;
    param.nu = 4;
    param.tol = tol;
    param.maxit = 1000;
    p2 = prox_l1(s,1,param);
    p3 = prox_tv1d(s,1,param); 
    
    if norm(p3-p2)/norm(p3)<1000*tol
        fprintf('  Test prox_l1 tv comp 1 OK\n')
    else
        fprintf('  Test prox_l1 tv comp Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end


function [errors]=test_fax(tol)

    errors = 0;
    
    N =100;
    
    x0 = 100*rand(N,1);
    
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
    
    param.verbose = 0;
    param.A = A;
    param.At = At;
    param.tight = 0;
    param.nu = 2*B;
    param.tol = tol;
    p2 = prox_l1(x0,1,param);
    
    paraml1.verbose = 0;
    f.prox = @(x,T) prox_l1(x,T,paraml1);
    f.eval = @(x) norm(x,1);
    param.f = f;
    
    p3 = prox_fax(x0,1,param); 
    
    if norm(p3-p2)/norm(p3)<1e-6
        fprintf('  Test prox_fac 1 OK\n');
    else
        fprintf('  Test prox_fac 1 Pas OK!!!!!!!!!!!!!!!!\n');
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end



function [errors]=test_l1_y()

    errors = 0;
    
    N =100;
    
    x0 = 100*rand(N,1);
    y = 100*rand(N,1);
    param.verbose = 0;
    param.y = y;
    param.tol = 1e-10;
    p2 = prox_l1(x0,1,param);
    
    
    p3 = prox_l1_admm(x0,1,param); 
    
    if norm(p3-p2)/norm(p3)<1e-8
        fprintf('  Test prox l1 y 1 OK\n');
    else
        fprintf('  Test prox l1 y 1 Pas OK!!!!!!!!!!!!!!!!\n');
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
 
  
    
    A = @(x) 1/sqrt(N)*fft(x);%  @(x) gradient_op1d(x);
    At = @(x) sqrt(N)*ifft(x);% @(x) -div_op1d(x);
    
    param.A = A;
    param.At = At;
    param.tight = 1;

    y = 100*rand(N,1);
    param.verbose = 0;
    param.y = y;
    param.tol = 1e-10;
    p2 = prox_l1(x0,1,param);
       
    p3 = prox_l1_admm(x0,1,param); 
    
    if norm(p3-p2)/norm(p3)<1e-8
        fprintf('  Test prox l1 y 2 OK\n');
    else
        fprintf('  Test prox l1 y 2 Pas OK!!!!!!!!!!!!!!!!\n');
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
    
       
    A = @(x) 1/sqrt(N)*fft(x);%  @(x) gradient_op1d(x);
    At = @(x) sqrt(N)*ifft(x);% @(x) -div_op1d(x);
    
    param.A = A;
    param.At = At;
    param.tight = 1;

    y = 100*rand(N,1);
    param.verbose = 0;
    param.y = y;
    param.tol = 1e-10;
    p2 = prox_l1(x0,1,param);
    param.tight = 0;
    p3 = prox_l1(x0,1,param); 
    
    if norm(p3-p2)/norm(p3)<1e-8
        fprintf('  Test prox l1 y 3 OK\n');
    else
        fprintf('  Test prox l1 y 3 Pas OK!!!!!!!!!!!!!!!!\n');
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
    
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
    
    param.verbose = 0;
    param.A = A;
    param.At = At;
    param.tight = 0;
    param.nu = 2*B;

    y = 100*rand(2*N,1);
    param.verbose = 0;
    param.y = y;
    param.tol = eps(10);
    p2 = prox_l1(x0,1,param);
    
    
    p3 = prox_l1_admm(x0,1,param); 
    
    if norm(p3-p2)/norm(p3)<1e-8
        fprintf('  Test prox l1 y 4 OK\n');
    else
        fprintf('  Test prox l1 y 4 Pas OK!!!!!!!!!!!!!!!!\n');
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end
