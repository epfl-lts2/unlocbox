function [ errors ] = test_tv( )
%TEST_TV This function test all the function related to tv
errors=0;

gsp_reset_seed(1)


errors=errors+test_norm1d(eps(10));
errors=errors+test_norm2d(eps(100));
errors=errors+test_norm3d(eps(100));
errors=errors+test_norm4d(eps(100));

errors=errors+test_prox_tv1d(eps(100));
errors=errors+test_prox_tv2d(eps(100));
errors=errors+test_prox_tv3d(eps(100));
errors=errors+test_prox_tv4d(eps(100));

errors=errors+test_prox_tv2d_mixed();

errors=errors+test_prox_tv2d_weights();
errors=errors+test_prox_tv3d_weights();
errors=errors+test_prox_tv4d_weights();

errors=errors+test_prox_tv1d_fast();





end

function [errors]=test_norm1d(tol)
    
    errors=0;

    s=rand(2,1);
    p2 = norm_tv1d(s);
    p3 = abs(s(1)-s(2)); 
    
    if norm(p3-p2)/norm(p3)<tol
        fprintf('  Test TV1d 1 OK\n')
    else
        fprintf('  Test TV1d 1 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
    s=rand(100,1);
    p2 = norm_tv1d(s);
    p3 = norm_tvnd(s); 
    
    if norm(p3-p2)/norm(p3)<tol
        fprintf('  Test TV1d 2 OK\n')
    else
        fprintf('  Test TV1d 2 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
    s=rand(100,2);
    p2 = norm_tv1d(s);
    p31 = norm_tv1d(s(:,1)); 
    p32 = norm_tv1d(s(:,2)); 
    p3 = [p31,p32]';
    
    if norm(p3-p2)/norm(p3)<tol
        fprintf('  Test TV1d 3 OK\n')
    else
        fprintf('  Test TV1d 3 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end


function [errors]=test_norm2d(tol)
    
    errors=0;

    s=zeros(10);
    s(45)=1;
    p2 = norm_tv(s);
    p3 = 2 + sqrt(2); 
    
    if norm(p3-p2)/norm(p3)<tol
        fprintf('  Test TV2d 1 OK\n')
    else
        fprintf('  Test TV2d 1 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
    s=rand(100);
    p2 = norm_tv(s);
    p3 = norm_tvnd(s); 
    
    if norm(p3-p2)/norm(p3)<tol
        fprintf('  Test TV2d 2 OK\n')
    else
        fprintf('  Test TV2d 2 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
    s=rand(100,100, 2);
    p2 = norm_tv(s);
    p31 = norm_tv(s(:,:,1)); 
    p32 = norm_tv(s(:,:,2)); 
    p3 = [p31,p32]';
    
    if norm(p3-p2)/norm(p3)<tol
        fprintf('  Test TV2d 3 OK\n')
    else
        fprintf('  Test TV2d 3 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end


function [errors]=test_norm3d(tol)
    
    errors=0;

    s=zeros(5,5,5);
    s(63)=1;
    p2 = norm_tv3d(s);
    p3 = 3 + sqrt(3); 
    
    if norm(p3-p2)/norm(p3)<tol
        fprintf('  Test TV3d 1 OK\n')
    else
        fprintf('  Test TV3d 1 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
    s=rand(10,10,10);
    p2 = norm_tv3d(s);
    p3 = norm_tvnd(s); 
    
    if norm(p3-p2)/norm(p3)<tol
        fprintf('  Test TV3d 2 OK\n')
    else
        fprintf('  Test TV3d 2 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
    s=rand(10,10,10, 2);
    p2 = norm_tv3d(s);
    p31 = norm_tv3d(s(:,:,:,1)); 
    p32 = norm_tv3d(s(:,:,:,2)); 
    p3 = [p31,p32]';
    
    if norm(p3-p2)/norm(p3)<tol
        fprintf('  Test TV3d 3 OK\n')
    else
        fprintf('  Test TV3d 3 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
end



function [errors]=test_norm4d(tol)
    
    errors=0;

    
    s=rand(10,10,10,10);
    p2 = norm_tv4d(s);
    p3 = norm_tvnd(s); 
    
    if norm(p3-p2)/norm(p3)<tol
        fprintf('  Test TV4d 1 OK\n')
    else
        fprintf('  Test TV4d 1 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
    s=rand(5,5,5,5, 2);
    p2 = norm_tv4d(s);
    p31 = norm_tv4d(s(:,:,:,:,1)); 
    p32 = norm_tv4d(s(:,:,:,:,2)); 
    p3 = [p31,p32]';
    
    if norm(p3-p2)/norm(p3)<tol
        fprintf('  Test TV4d 2 OK\n')
    else
        fprintf('  Test TV4d 2 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
end



function [errors]=test_prox_tv1d(tol)
    
    errors=0;
    
    N = 100;
    
    so = zeros(N,1);
    so(1:50) = 1;
    
    nstv = norm_tv1d(so);
    
    noise = 0.05*(rand(N,1)-0.5);
    noise(1:50) = noise(1:50)-mean(noise(1:50));
    noise(51:100) = noise(51:100)-mean(noise(51:100));
    
    nn2 = norm(noise)^2;
    
    lambda = nn2/nstv*20;
    
    s= so + noise;
    
    param.verbose = 0;
    param.tol = tol;
    p2 = prox_tv1d(s,lambda,param);
    p3 = so; 
    
    if abs(norm_tv1d(p3)-norm_tv1d(p2))<0.05
        fprintf('  Test prox_tv1d 1 OK\n')
    else
        fprintf('  Test prox_tv1d 1 Pas OK!!!!!!!!!!!!!!!!\n')
        errors= errors +1;
        plot(1:N,so,1:N,p2,'r',1:N,s,'g')
    end
    
    
    s=rand(100, 2);
    p2 = prox_tv1d(s,lambda,param);
    p31 = prox_tv1d(s(:,1),lambda,param); 
    p32 = prox_tv1d(s(:,2),lambda,param); 
    p3 = [p31,p32];
    
     if norm(p3(:)-p2(:))/norm(p3(:))<100*tol
        fprintf('  Test prox_tv1d 2 OK\n')
    else
        fprintf('  Test prox_tv1d 2 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3(:)-p2(:))/norm(p3)
        errors= errors +1;
     end
    
    s=ones(5, 2);
    p2 = prox_tv1d(s,1,param);
    p3 = s;
    
     if norm(p3-p2)/norm(p3)<tol
        fprintf('  Test prox_tv1d 3 OK\n')
    else
        fprintf('  Test prox_tv1d 3 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end



function [errors]=test_prox_tv2d(tol)
    
    errors=0;
    
    N = 10;
    
    so = zeros(N,N);
    so(1:50) = 1;
    
    nstv = norm_tv(so);
    
    noise = 0.05*(rand(N,N)-0.5);
    noise(1:50) = noise(1:50)-mean(noise(1:50));
    noise(51:100) = noise(51:100)-mean(noise(51:100));
    
    nn2 = norm(noise(:))^2;
    
    lambda = nn2/nstv*10;
    
    s= so + noise;
    
    param.verbose = 0;
    param.tol = tol;
    p2 = prox_tv(s,lambda,param);
    p3 = so; 
    
    if abs(norm_tv(p3)-norm_tv(p2))<0.1
        fprintf('  Test prox_tv2d 1 OK\n')
    else
        fprintf('  Test prox_tv2d 1 Pas OK!!!!!!!!!!!!!!!!\n')
        errors= errors +1;
        figure
        imagesc(so)
        figure
        imagesc(p3)
        figure
      imagesc(s)

    end
    
    
     s=rand(5,5, 2);
    p2 = prox_tv(s,lambda,param);
    p31 = prox_tv(s(:,:,1),lambda,param); 
    p32 = prox_tv(s(:,:,2),lambda,param); 
    p3 = cat(3,p31,p32);
    
     if norm(p3(:)-p2(:))/norm(p3(:))<1e-7
        fprintf('  Test prox_tv2d 2 OK\n')
    else
        fprintf('  Test prox_tv2d 2 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3(:)-p2(:))/norm(p3(:))
        errors= errors +1;
     end
    
    s=ones(10,10,3);
    p2 = prox_tv(s,1,param);
    p3 = s;
    
     if norm(p3(:)-p2(:))/norm(p3(:))<tol
        fprintf('  Test prox_tv2d 3 OK\n')
    else
        fprintf('  Test prox_tv2d 3 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3(:)-p2(:))/norm(p3(:))
        errors= errors +1;
    end
end


function [errors]=test_prox_tv3d(tol)
    
    errors=0;
    
    N = 5;
    
    so = zeros(N,N,N);
    so(1:50) = 1;
    
    nstv = norm_tv3d(so);
    
    noise = 0.05*(rand(N,N,N)-0.5);
    noise(1:50) = noise(1:50)-mean(noise(1:50));
    noise(51:125) = noise(51:125)-mean(noise(51:125));
    
    nn2 = norm(noise(:))^2;
    
    lambda = nn2/nstv*10;
    
    s= so + noise;
    
    param.verbose = 0;
    param.tol = tol;
    p2 = prox_tv3d(s,lambda,param);
    p3 = so; 
    
    if abs(norm_tv3d(p3)-norm_tv3d(p2))/norm_tv3d(p2)<0.01
        fprintf('  Test prox_tv3d 1 OK\n')
    else
        fprintf('  Test prox_tv3d 1 Pas OK!!!!!!!!!!!!!!!!\n')
        errors= errors +1;
        abs(norm_tv3d(p3)-norm_tv3d(p2))/norm_tv3d(p2)

    end
    
    
     s=rand(5,5,5, 2);
    p2 = prox_tv3d(s,lambda,param);
    p31 = prox_tv3d(s(:,:,:,1),lambda,param); 
    p32 = prox_tv3d(s(:,:,:,2),lambda,param); 
    p3 = cat(4,p31,p32);
    
     if norm(p3(:)-p2(:))/norm(p3(:))<1e-7
        fprintf('  Test prox_tv3d 2 OK\n')
    else
        fprintf('  Test prox_tv3d 2 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3(:)-p2(:))/norm(p3(:))
        errors= errors +1;
     end
    
    s=ones(5,5,5,3);
    p2 = prox_tv3d(s,1,param);
    p3 = s;
    
     if norm(p3(:)-p2(:))/norm(p3(:))<tol
        fprintf('  Test prox_tv3d 3 OK\n')
    else
        fprintf('  Test prox_tv3d 3 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3(:)-p2(:))/norm(p3(:))
        errors= errors +1;
    end
end




function [errors]=test_prox_tv4d(tol)
    
    errors=0;
    
    N = 5;
    
    so = zeros(N,N,N,N);
    so(1:300) = 1;
    
    nstv = norm_tv4d(so);
    
    noise = 0.05*(rand(N,N,N,N)-0.5);
    noise(1:300) = noise(1:300)-mean(noise(1:300));
    noise(301:625) = noise(301:625)-mean(noise(301:625));
    
    nn2 = norm(noise(:))^2;
    
    lambda = nn2/nstv*10;
    
    s= so + noise;
    
    param.verbose = 0;
    param.tol = tol;
    p2 = prox_tv4d(s,lambda,param);
    p3 = so; 
    
    if abs(norm_tv4d(p3)-norm_tv4d(p2))/norm_tv4d(p2)<0.01
        fprintf('  Test prox_tv4d 1 OK\n')
    else
        fprintf('  Test prox_tv4d 1 Pas OK!!!!!!!!!!!!!!!!\n')
        errors= errors +1;
        abs(norm_tv4d(p3)-norm_tv4d(p2))/norm_tv4d(p2)

    end
    
    
     s=rand(5,5,5,5, 2);
    p2 = prox_tv4d(s,lambda,param);
    p31 = prox_tv4d(s(:,:,:,:,1),lambda,param); 
    p32 = prox_tv4d(s(:,:,:,:,2),lambda,param); 
    p3 = cat(5,p31,p32);
    
     if norm(p3(:)-p2(:))/norm(p3(:))<1e-7
        fprintf('  Test prox_tv4d 2 OK\n')
    else
        fprintf('  Test prox_tv4d 2 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3(:)-p2(:))/norm(p3(:))
        errors= errors +1;
     end
    
    s=ones(5,5,5,5,3);
    p2 = prox_tv4d(s,1,param);
    p3 = s;
    
     if norm(p3(:)-p2(:))/norm(p3(:))<tol
        fprintf('  Test prox_tv4d 3 OK\n')
    else
        fprintf('  Test prox_tv4d 3 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3(:)-p2(:))/norm(p3(:))
        errors= errors +1;
    end
end



function [errors]=test_prox_tv2d_mixed()

    errors = 0;
    
    N = 10;
    s = 10*rand(N,N);

    
    param.verbose = 0;
    param.A = @(x) grad2d_help(x);
    param.At = @(x) div2d_help(x);
    param.tight = 0;
    param.nu = 4;
    param.maxit = 1000;
    param.tol = eps(100);
    g_d = zeros(1,2*N^2);
    for ii = 1:N^2
        g_d(2*ii-1) = ii;
        g_d(2*ii) = ii+N^2;
    end
    paraml21.verbose = 0;
    paraml21.g_d = g_d;
    paraml21.g_t = 2*ones(1,N^2);
    
    f.prox = @(x,T) prox_l21(x,T,paraml21);
    f.eval = @(x)   norm_l21(x,paraml21.g_d,paraml21.g_t);
    param.f = f;
    
    p2 = prox_fax(s(:),1,param);
    p3 = prox_tv(s,1,param); 
    
    if norm(p3(:)-p2(:))/norm(p3(:))<1e-6
        fprintf('  Test prox_l1 tv comp 1 OK\n')
    else
        fprintf('  Test prox_l1 tv comp Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3(:)-p2(:))/norm(p3(:))
        errors= errors +1;
    end
    
end



function [errors]=test_prox_tv2d_weights()
    
    errors=0;
    
    N = 10;
    
    s = rand(N,N);
   
    lambda = 0.1;
    
    param.verbose = 0;
    param.tol = eps;
    param.weights = [1,0];
    p2 = prox_tv(s,lambda,param);
    p3 = prox_tv1d(s,lambda,param);
    
    if norm(p3(:)-p2(:))/norm(p3(:))<1e-7
        fprintf('  Test prox_tv2d weights 1 OK\n')
    else
        fprintf('  Test prox_tv2d weights 1 Pas OK!!!!!!!!!!!!!!!!\n')
        errors= errors +1;
        norm(p3(:)-p2(:))/norm(p3(:))
    end
    


    
    N = 10;
    s = 10*rand(N,N);

    
    wx = 0.3;
    wy = 0.8;
    
   
    
    param.verbose = 0;
    param.A = @(x) grad2d_help(x,wx,wy);
    param.At = @(x) div2d_help(x,wx,wy);
    param.tight = 0;
    param.nu = 4;
    param.maxit = 1000;
    param.tol = eps(100);
    
    
    g_d = zeros(1,2*N^2);
    
    for ii = 1:N^2
        g_d(2*ii-1) = ii;
        g_d(2*ii) = ii+N^2;
    end
    paraml21.verbose = 0;
    paraml21.g_d = g_d;
    paraml21.g_t = 2*ones(1,N^2);
    
    f.prox = @(x,T) prox_l21(x,T,paraml21);
    f.eval = @(x)   norm_l21(x,paraml21.g_d,paraml21.g_t);
    param.f = f;
    
    p2 = prox_fax(s(:),1,param);
    param.weights = [wx,wy];
    p3 = prox_tv(s,1,param); 
    
    if norm(p3(:)-p2(:))/norm(p3(:))<1e-5
        fprintf('  Test prox_tv2d weight 2 OK\n')
    else
        fprintf('  Test prox_tv2d weight 2 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3(:)-p2(:))/norm(p3(:))
        errors= errors +1;
    end
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    N = 10;
    s = 10*rand(N,N);

    
    wx = 0.9;
    wy = 5;
    
   
    
    param.verbose = 0;
    param.A = @(x) grad2d_help(x,wx,wy);
    param.At = @(x) div2d_help(x,wx,wy);
    param.tight = 0;
    param.nu = 4*wy;
    param.maxit = 1000;
    param.tol = eps(100);
    
    
    g_d = zeros(1,2*N^2);
    
    for ii = 1:N^2
        g_d(2*ii-1) = ii;
        g_d(2*ii) = ii+N^2;
    end
    paraml21.verbose = 0;
    paraml21.g_d = g_d;
    paraml21.g_t = 2*ones(1,N^2);
    
    f.prox = @(x,T) prox_l21(x,T,paraml21);
    f.eval = @(x)   norm_l21(x,paraml21.g_d,paraml21.g_t);
    param.f = f;
    
    p2 = prox_fax(s(:),1,param);
    param.weights = [wx,wy];
    p3 = prox_tv(s,1,param); 
    
    if norm(p3(:)-p2(:))/norm(p3(:))<4e-4
        fprintf('  Test prox_tv2d weight 3 OK\n')
    else
        fprintf('  Test prox_tv2d weight 3 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3(:)-p2(:))/norm(p3(:))
        errors= errors +1;
    end
    
    
  
end




function [errors]=test_prox_tv3d_weights()
    
    errors=0;
    
    N = 5;
    
    s = rand(N,N,N);
   
    lambda = 0.1;
    
    param.verbose = 0;
    param.tol = eps;  
    param.maxit = 1000;
    p2 = prox_tv(s,lambda,param);
    param.weights = [1,1,0];
    p3 = prox_tv3d(s,lambda,param);
    
    if norm(p3(:)-p2(:))/norm(p3(:))<1e-7
        fprintf('  Test prox_tv3d weights 1 OK\n')
    else
        fprintf('  Test prox_tv3d weights 1 Pas OK!!!!!!!!!!!!!!!!\n')
        errors= errors +1;
        norm(p3(:)-p2(:))/norm(p3(:))
    end
    


    
    N = 5;
    s = 10*rand(N,N,N);

    
    wx = 0.3;
    wy = 0.8;
    wz = 0.5;
    
   
    
    param.verbose = 0;
    param.A = @(x) grad3d_help(x,wx,wy,wz);
    param.At = @(x) div3d_help(x,wx,wy,wz);
    param.tight = 0;
    param.nu = 4;
    param.maxit = 1000;
    param.tol = eps(100);
    
    
    g_d = zeros(1,3*N^3);
    
    for ii = 1:N^3
        g_d(3*ii-2) = ii;
        g_d(3*ii-1) = ii+N^3;
        g_d(3*ii) = ii+2*N^3;
    end
    paraml21.verbose = 0;
    paraml21.g_d = g_d;
    paraml21.g_t = 3*ones(1,N^3);
    
    f.prox = @(x,T) prox_l21(x,T,paraml21);
    f.eval = @(x)   norm_l21(x,paraml21.g_d,paraml21.g_t);
    param.f = f;
    
    p2 = prox_fax(s(:),1,param);
    param.weights = [wx,wy,wz];
    p3 = prox_tv3d(s,1,param); 
    
    if norm(p3(:)-p2(:))/norm(p3(:))<1e-6
        fprintf('  Test prox_tv3d weight 2 OK\n')
    else
        fprintf('  Test prox_tv3d weight 2 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3(:)-p2(:))/norm(p3(:))
        errors= errors +1;
    end
  
end



function [errors]=test_prox_tv4d_weights()
    
    errors=0;
    
    N = 5;
    
    s = rand(N,N,N,N);
   
    lambda = 0.1;
    
    param.verbose = 0;
    param.tol = eps;  
    param.maxit = 700;
    p2 = prox_tv3d(s,lambda,param);
    param.weights = [1,1,1,0];
    p3 = prox_tv4d(s,lambda,param);
    
    if norm(p3(:)-p2(:))/norm(p3(:))<1e-5
        fprintf('  Test prox_tv4d weights 1 OK\n')
    else
        fprintf('  Test prox_tv4d weights 1 Pas OK!!!!!!!!!!!!!!!!\n')
        errors= errors +1;
        norm(p3(:)-p2(:))/norm(p3(:))
    end
    


    
    N = 4;
    s = 5+rand(N,N,N,N);

    
    wx = 0.3;
    wy = 0.8;
    wz = 0.5;
    wt = 0.6;

    
   
    
    param.verbose = 0;
    param.A = @(x) grad4d_help(x,wx,wy,wz,wt);
    param.At = @(x) div4d_help(x,wx,wy,wz,wt);
    param.tight = 0;
    param.nu = 4;
    param.maxit = 2000;
    param.tol = eps(100);
    
    
    g_d = zeros(1,3*N^4);
    
    for ii = 1:N^4
        g_d(4*ii-3) = ii;
        g_d(4*ii-2) = ii+N^4;
        g_d(4*ii-1) = ii+2*N^4;
        g_d(4*ii)   = ii+3*N^4;
    end
    paraml21.verbose = 0;
    paraml21.g_d = g_d;
    paraml21.g_t = 4*ones(1,N^4);
    
    f.prox = @(x,T) prox_l21(x,T,paraml21);
    f.eval = @(x)   norm_l21(x,paraml21.g_d,paraml21.g_t);
    param.f = f;
    
    p2 = prox_fax(s(:),1,param);
    param.weights = [wx,wy,wz,wt];
    p3 = prox_tv3d(s,1,param); 
    
    if norm(p3(:)-p2(:))/norm(p3(:))<3e-2
        fprintf('  Test prox_tv4d weight 2 OK\n')
    else
        fprintf('  Test prox_tv4d weight 2 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3(:)-p2(:))/norm(p3(:))
        errors= errors +1;
    end
  
end




function g = grad2d_help(x,wx,wy)

    N = numel(x);
    x = reshape(x,sqrt(N),sqrt(N));
    if nargin > 1
        [a,b]=gradient_op(x,wx,wy);
    else
        [a,b]=gradient_op(x);
    end
    g = reshape(cat(3,a,b),1,[]);
end

function g = div2d_help(x,wx,wy)
    N = length(x);
    x = reshape(x,sqrt(N/2),sqrt(N/2),2);
    if nargin > 1
        g = -reshape(div_op(x(:,:,1),x(:,:,2),wx,wy),[],1);
    else
        g = -reshape(div_op(x(:,:,1),x(:,:,2)),[],1);
    end
end




function g = grad3d_help(x,wx,wy,wz)

    N = numel(x);
    x = reshape(x,round(N^(1/3)),round(N^(1/3)),round(N^(1/3)));
    if nargin > 1
        [a,b,c]=gradient_op3d(x,wx,wy,wz);
    else
        [a,b,c]=gradient_op3d(x);
    end
    g = reshape(cat(4,a,b,c),1,[]);
end

function g = div3d_help(x,wx,wy,wz)
    N = numel(x);
    x = reshape(x,round((N/3)^(1/3)),round((N/3)^(1/3)),round((N/3)^(1/3)),3);
    if nargin > 1
        g = -reshape(div_op3d(x(:,:,:,1),x(:,:,:,2),x(:,:,:,3),wx,wy,wz),[],1);
    else
        g = -reshape(div_op3d(x(:,:,:,1),x(:,:,:,2),x(:,:,:,3)),[],1);
    end
end




function g = grad4d_help(x,wx,wy,wz,wt)

    N = numel(x);
    x = reshape(x,round(N^(1/4)),round(N^(1/4)),round(N^(1/4)),round(N^(1/4)));
    if nargin > 1
        [a,b,c,d]=gradient_op4d(x,wx,wy,wz,wt);
    else
        [a,b,c,d]=gradient_op4d(x);
    end
    g = reshape(cat(5,a,b,c,d),1,[]);
end

function g = div4d_help(x,wx,wy,wz,wt)
    N = numel(x);
    x = reshape(x,round((N/4)^(1/4)),round((N/4)^(1/4)),round((N/4)^(1/4)),round((N/4)^(1/4)),4);
    if nargin > 1
        g = -reshape(div_op4d(x(:,:,:,:,1),x(:,:,:,:,2),x(:,:,:,:,3),x(:,:,:,:,4),wx,wy,wz,wt),[],1);
    else
        g = -reshape(div_op4d(x(:,:,:,:,1),x(:,:,:,:,2),x(:,:,:,:,3),x(:,:,:,:,4)),[],1);
    end
end





function [errors]=test_prox_tv1d_fast()
    
    errors=0;
    
    N = 100;
    
    s = rand(N,1);
    lambda = 0.1;
    param.use_fast = 0;
    param.tol = 1e-16;
    param.maxit = 1000;
    param.verbose = 0;
    p2 = prox_tv1d(s,lambda,param);
    param.use_fast = 1;
    p3 = prox_tv1d(s,lambda,param);
   
    
     if norm(p3(:)-p2(:))/norm(p3(:))<1e-8
        fprintf('  Test prox_tv1d fast OK\n')
    else
        fprintf('  Test prox_tv1d fast Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3(:)-p2(:))/norm(p3)
        errors= errors +1;
     end
    
     
    s = rand(N,10);
    lambda = 0.1;
    param.use_fast = 0;
    param.tol = 1e-16;
    param.maxit = 1000;
    p2 = prox_tv1d(s,lambda,param);
    param.use_fast = 1;
%     for ii = 1:10
%         p3(:,ii) = prox_tv1d(s(:,ii),lambda,param);
%     end
   p3 = prox_tv1d(s,lambda,param);
    
     if norm(p3(:)-p2(:))/norm(p3(:))<1e-8
        fprintf('  Test prox_tv1d fast Ndim 2OK\n')
    else
        fprintf('  Test prox_tv1d fast Ndim Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3(:)-p2(:))/norm(p3)
        errors= errors +1;
     end
    
end
