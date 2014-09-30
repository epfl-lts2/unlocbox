function [b]=test_mixed_sparsity()
%% Test the mixed norm
epsilon=10e-12;
b=0;

b=b+test_soft(epsilon);
b=b+test_softB(epsilon);

b=b+test_l21norm(epsilon);
b=b+test_l21norm_simple(epsilon);
b=b+test_l21prox(epsilon);
b=b+test_l21norm_matrix(epsilon);
b=b+test_l21prox_matrix(epsilon);

b=b+test_find_mg(epsilon);

b=b+test_l12norm(epsilon);
b=b+test_l12norm_simple(epsilon);
b=b+test_l12norm_matrix(epsilon);
b=b+test_l12prox(epsilon);
b=b+test_l12prox_matrix(epsilon);

b=b+test_linf1norm(epsilon);
b=b+test_linf1norm_simple(epsilon);
b=b+test_linf1prox(epsilon);
b=b+test_linf1norm_matrix(epsilon);
b=b+test_linf1prox_matrix(epsilon);

end

%%%%%%%%%%%%%%%%%%%%%%
%  -- test soft     --
%%%%%%%%%%%%%%%%%%%%%%
function [b] = test_soft(epsilon)
    L=100;    
    T=1;
    x=5*randn(L,1);
    soft= @(z,T) sign(z).*max(abs(z)-T, 0);
    if (norm(soft_threshold(x,T)-soft(x,T))<epsilon)
        fprintf('Test Soft threshold - ok\n')
        b=0;
    else
        fprintf('Test Soft threshold - Error !!!!!!!!!!!!!!\n')
        b=1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%
%  -- test softB    --
%%%%%%%%%%%%%%%%%%%%%%
function [b] = test_softB(epsilon)
    L=100;    
    T=1;
    x=5*randn(L,1);
    softB = @(z, T) sign(z).*max(abs(z)-T.*abs(z), 0);   
    if (norm(soft_thresholdb(x,T)-softB(x,T))<epsilon)
        fprintf('Test Soft threshold b - ok\n')
        b=0;
    else
        fprintf('Test Soft threshold b - Error !!!!!!!!!!!!!!\n')
        b=1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%
%  -- test norm L21 --
%%%%%%%%%%%%%%%%%%%%%%
function [b] = test_l21norm(epsilon)
    x= rand(6,1);
    p=randperm(6);
    param.g_d=p;
    param.g_t=[2 2 2];
    x1=[x(p(1)),x(p(2))];
    x2=[x(p(3)),x(p(4))];
    x3=[x(p(5)),x(p(6))];
    nx=norm([norm(x1),norm(x2),norm(x3)],1);
    %norm_l21(x,param.g_d,param.g_t)
    if (norm(nx-norm_l21(x,param.g_d,param.g_t))<epsilon)
        fprintf('Test L21 norm - Group same size - ok\n')
        b=0;
    else
        fprintf('Test L21 norm - Group same size - Error !!!!!!!!!!!!!!\n')
        b=1;
    end
     
    x= rand(6,1);
    p=randperm(6);
    param.g_d=p;
    param.g_t=[2 4];
    x1=[x(p(1)),x(p(2))];
    x2=[x(p(3)),x(p(4)),x(p(5)),x(p(6))];
    nx=norm([norm(x1),norm(x2)],1);
    % norm_l21(x,param.g_d,param.g_t);
    if (norm(nx-norm_l21(x,param.g_d,param.g_t))<epsilon)
        fprintf('Test L21 norm - Group diff size - ok\n')
        b=0+b;
    else
        fprintf('Test L21 norm - Group diff size - Error !!!!!!!!!!!!!!\n')
        b=1+b;
    end
end


%%%%%%%%%%%%%%%%%%%%%%
%  -- test norm L21 simple --
%%%%%%%%%%%%%%%%%%%%%%
function [b] = test_l21norm_simple(epsilon)
    x= rand(1,100);
    nx = norm(x);
    if abs(nx-norm_l21(x))<epsilon
        fprintf('Test L21 norm - simple 1 - ok\n')
        b=0;
    else
        fprintf('Test L21 norm - simple 1 - Error !!!!!!!!!!!!!!\n')
        b=1;
    end
     
    x= rand(100,1);
    nx = norm(x,1);
    if abs(nx-norm_l21(x))<epsilon
        fprintf('Test L21 norm - simple 2 - ok\n')
        b=0+b;
    else
        fprintf('Test L21 norm - simple 2 - Error !!!!!!!!!!!!!!\n')
        b=1+b;
    end
end


%%%%%%%%%%%%%%%%%%%%%%
%  -- test prox L21 --
%%%%%%%%%%%%%%%%%%%%%%
function [b] = test_l21prox(epsilon)

    softB = @(z, T) sign(z).*max(abs(z)-T.*abs(z), 0);   
    lambda=0.5;
    x= rand(6,1);
    p=randperm(6);
    param.g_d=p;
    param.g_t=[2 2 2];
    x1=[x(p(1)),x(p(2))];
    x2=[x(p(3)),x(p(4))];
    x3=[x(p(5)),x(p(6))];
    param.verbose=0;
    n1=softB(x1,lambda/norm(x1,2));
    n2=softB(x2,lambda/norm(x2,2));
    n3=softB(x3,lambda/norm(x3,2));
    n=[n1,n2,n3];
    nx=n(back_perm(p'))';
    %prox_l21(x,lambda,param);
    if (norm(nx- prox_l21(x,lambda,param))<epsilon)
        fprintf('Test L21 prox - Group same size - ok\n')
        b=0;
    else
        fprintf('Test L21 prox - Group same size - Error !!!!!!!!!!!!!!\n')
        b=1;
    end
     
    x= rand(6,1);
    p=randperm(6);

    param.g_d=p;
    param.g_t=[2 4];
    x1=[x(p(1)),x(p(2))];
    x2=[x(p(3)),x(p(4)),x(p(5)),x(p(6))];
    n1=softB(x1,lambda/norm(x1,2));
    n2=softB(x2,lambda/norm(x2,2));
    n=[n1,n2];
    nx=n(back_perm(p'))';
    if (norm(nx- prox_l21(x,lambda,param))<epsilon)
        fprintf('Test L21 prox - Group diff size - ok\n')
        b=0+b;
    else
        fprintf('Test L21 prox - Group diff size - Error !!!!!!!!!!!!!!\n')
        b=1+b;
    end
 
end


%%%%%%%%%%%%%%%%%%%%%%
%  -- test norm Linf1 --
%%%%%%%%%%%%%%%%%%%%%%
function [b] = test_linf1norm(epsilon)
    x= rand(6,1);
    p=randperm(6);
    param.g_d=p;
    param.g_t=[2 2 2];
    x1=[x(p(1)),x(p(2))];
    x2=[x(p(3)),x(p(4))];
    x3=[x(p(5)),x(p(6))];
    nx=norm([norm(x1,Inf),norm(x2,Inf),norm(x3,Inf)],1);

    if (norm(nx-norm_linf1(x,param.g_d,param.g_t))<epsilon)
        fprintf('Test Linf1 norm - Group same size - ok\n')
        b=0;
    else
        fprintf('Test Linf1 norm - Group same size - Error !!!!!!!!!!!!!!\n')
        b=1;
    end
     
    x= rand(6,1);
    p=randperm(6);
    param.g_d=p;
    param.g_t=[2 4];
    x1=[x(p(1)),x(p(2))];
    x2=[x(p(3)),x(p(4)),x(p(5)),x(p(6))];
    nx=norm([norm(x1,Inf),norm(x2,Inf)],1);

    if (norm(nx-norm_linf1(x,param.g_d,param.g_t))<epsilon)
        fprintf('Test Linf1 norm - Group diff size - ok\n')
        b=0+b;
    else
        fprintf('Test Linf1 norm - Group diff size - Error !!!!!!!!!!!!!!\n')
        b=1+b;
    end
end


%%%%%%%%%%%%%%%%%%%%%%
%  -- test norm Linf1 simple --
%%%%%%%%%%%%%%%%%%%%%%
function [b] = test_linf1norm_simple(epsilon)
    x= rand(1,100);
    nx = max(abs(x));
    if abs(nx-norm_linf1(x))<epsilon
        fprintf('Test Linf1 norm - simple 1 - ok\n')
        b=0;
    else
        fprintf('Test Linf1 norm - simple 1 - Error !!!!!!!!!!!!!!\n')
        b=1;
    end
     
    x= rand(100,1);
    nx = norm(x,1);
    if abs(nx-norm_linf1(x))<epsilon
        fprintf('Test Linf1 norm - simple 2 - ok\n')
        b=0+b;
    else
        fprintf('Test Linf1 norm - simple 2 - Error !!!!!!!!!!!!!!\n')
        b=1+b;
    end
end

%%%%%%%%%%%%%%%%%%%%%%
%  -- test prox Linf1 --
%%%%%%%%%%%%%%%%%%%%%%
function [b] = test_linf1prox(epsilon)
    softB = @(z, T) sign(z).*max(abs(z)-T.*abs(z), 0);   
    
    lambda=0.5;
    x= rand(6,1);
    p=randperm(6);
    param.g_d=p;
    param.g_t=[2 2 2];
    x1=[x(p(1)),x(p(2))];
    x2=[x(p(3)),x(p(4))];
    x3=[x(p(5)),x(p(6))];
    param.verbose=0;
    n1=softB(x1,lambda/norm(x1,Inf));
    n2=softB(x2,lambda/norm(x2,Inf));
    n3=softB(x3,lambda/norm(x3,Inf));
    n=[n1,n2,n3];
    nx=n(back_perm(p'))';

    if (norm(nx- prox_linf1(x,lambda,param))<epsilon)
        fprintf('Test Linf1 prox - Group same size - ok\n')
        b=0;
    else
        fprintf('Test Linf1 prox - Group same size - Error !!!!!!!!!!!!!!\n')
        b=1;
    end
     
    x= rand(6,1);
    p=randperm(6);

    param.g_d=p;
    param.g_t=[2 4];
    x1=[x(p(1)),x(p(2))];
    x2=[x(p(3)),x(p(4)),x(p(5)),x(p(6))];
    n1=softB(x1,lambda/norm(x1,Inf));
    n2=softB(x2,lambda/norm(x2,Inf));
    n=[n1,n2];
    nx=n(back_perm(p'))';
    if (norm(nx- prox_linf1(x,lambda,param))<epsilon)
        fprintf('Test Linf1 prox - Group diff size - ok\n')
        b=0+b;
    else
        fprintf('Test Linf1 prox - Group diff size - Error !!!!!!!!!!!!!!\n')
        b=1+b;
    end
 
end



function [b] = test_l21norm_matrix(epsilon)

    x= rand(6,10);
    
    nx=sum(sqrt(sum(abs(x).^2,2)),1);
    
    

    if abs(nx-norm_l21(x))<epsilon
        fprintf('Test L21 norm - matrix - ok\n')
        b=0;
    else
        fprintf('Test L21 norm - matrix - Error !!!!!!!!!!!!!!\n')
        b=1;
    end


end

function [b] = test_linf1norm_matrix(epsilon)

    x= rand(6,10);  
    nx=sum(max(abs(x),[],2),1);
      
    if abs(nx-norm_linf1(x))<epsilon
        fprintf('Test Linf1 norm - matrix - ok\n')
        b=0;
    else
        fprintf('Test Linf1 norm - matrix - Error !!!!!!!!!!!!!!\n')
        b=1;
    end


end


function [b] = test_l21prox_matrix(epsilon)

    b = 0;

    softB = @(z, T) sign(z).*max(abs(z)-T.*abs(z), 0);   
    lambda=0.5;

  
%     n=softB(x,lambda./repmat(sqrt(sum(abs(x).^2,2)),1,size(x,2)));
    param.verbose = 0;
    
    x = rand(3,2);
    
    x1=x(1,:);
    x2=x(2,:);
    x3=x(3,:);
    param.verbose=0;
    n1=softB(x1,lambda/norm(x1,2));
    n2=softB(x2,lambda/norm(x2,2));
    n3=softB(x3,lambda/norm(x3,2));
    
    n = [n1 ; n2; n3];
    
    n2 = prox_l21(x,lambda,param);



    if norm(n(:) - n2(:)) <epsilon
        fprintf('Test L21 prox - matrix 1 - ok\n')
    else
        fprintf('Test L21 prox - matrix 1 - Error !!!!!!!!!!!!!!\n')
        b=b+1;
    end
    
    
    
    x = 5*rand(100,1);


    param.verbose = 0;

  
    
    n2 = prox_l21(x,lambda,param);
    n = prox_l1(x,lambda,param);


    if norm(n(:) - n2(:)) <epsilon
        fprintf('Test L21 prox - matrix 2 - ok\n')
    else
        fprintf('Test L21 prox - matrix 2 - Error !!!!!!!!!!!!!!\n')
        b=b+1;
    end
    
    

end




function [b] = test_linf1prox_matrix(epsilon)

    b = 0;
    softB = @(z, T) sign(z).*max(abs(z)-T.*abs(z), 0);   
    lambda=0.5;

  
%     n=softB(x,lambda./repmat(sqrt(sum(abs(x).^2,2)),1,size(x,2)));
    param.verbose = 0;
    
    x = rand(3,2);
    
    x1=x(1,:);
    x2=x(2,:);
    x3=x(3,:);
    param.verbose=0;
    n1=softB(x1,lambda/norm(x1,Inf));
    n2=softB(x2,lambda/norm(x2,Inf));
    n3=softB(x3,lambda/norm(x3,Inf));
    
    n = [n1 ; n2; n3];
    
    n2 = prox_linf1(x,lambda,param);



    if norm(n(:) - n2(:)) <epsilon
        fprintf('Test Linf1 prox - matrix 1 - ok\n')
    else
        fprintf('Test Linf1 prox - matrix 1 - Error !!!!!!!!!!!!!!\n')
        b=b+1;
    end
    
    
    
    x = 5*rand(100,1);


    param.verbose = 0;

  
    
    n2 = prox_linf1(x,lambda,param);
    n = prox_l1(x,lambda,param);


    if norm(n(:) - n2(:)) <epsilon
        fprintf('Test Linf1 prox - matrix 2 - ok\n')
    else
        fprintf('Test Linf1 prox - matrix 2 - Error !!!!!!!!!!!!!!\n')
        b=b+1;
    end

end






%%%%%%%%%%%%%%%%%%%%%%
%  -- test norm L12 --
%%%%%%%%%%%%%%%%%%%%%%
function [b] = test_l12norm(epsilon)
    x= rand(6,1);
    p=randperm(6);
    param.g_d=p;
    param.g_t=[2 2 2];
    x1=[x(p(1)),x(p(2))];
    x2=[x(p(3)),x(p(4))];
    x3=[x(p(5)),x(p(6))];
    nx=norm([norm(x1,1),norm(x2,1),norm(x3,1)],2);
    if (norm(nx-norm_l12(x,param.g_d,param.g_t))<epsilon)
        fprintf('Test L12norm - Group same size - ok\n')
        b=0;
    else
        fprintf('Test L12 norm - Group same size - Error !!!!!!!!!!!!!!\n')
        b=1;
    end
     
    x= rand(6,1);
    p=randperm(6);
    param.g_d=p;
    param.g_t=[2 4];
    x1=[x(p(1)),x(p(2))];
    x2=[x(p(3)),x(p(4)),x(p(5)),x(p(6))];
    nx=norm([norm(x1,1),norm(x2,1)],2);
    if (norm(nx-norm_l12(x,param.g_d,param.g_t))<epsilon)
        fprintf('Test L12 norm - Group diff size - ok\n')
        b=0+b;
    else
        fprintf('Test L12 norm - Group diff size - Error !!!!!!!!!!!!!!\n')
                norm(nx-norm_l12(x,param.g_d,param.g_t))
        b=1+b;
    end
end


%%%%%%%%%%%%%%%%%%%%%%
%  -- test norm L12 simple --
%%%%%%%%%%%%%%%%%%%%%%
function [b] = test_l12norm_simple(epsilon)
    x= rand(1,100);
    nx = norm(x,1);
    if abs(nx-norm_l12(x))<epsilon
        fprintf('Test L12 norm - simple 1 - ok\n')
        b=0;
    else
        fprintf('Test L12 norm - simple 1 - Error !!!!!!!!!!!!!!\n')
        b=1;
    end
     
    x= rand(100,1);
    nx = norm(x,2);
    if abs(nx-norm_l12(x))<epsilon
        fprintf('Test L12 norm - simple 2 - ok\n')
        b=0+b;
    else
        fprintf('Test L12 norm - simple 2 - Error !!!!!!!!!!!!!!\n')
        b=1+b;
    end
end




%%%%%%%%%%%%%%%%%%%%%%
%  -- test prox L12 --
%%%%%%%%%%%%%%%%%%%%%%
function [b] = test_l12prox(epsilon)

    softC = @(z, T) sign(z).*max(abs(z)-T, 0);   
    lambda=0.5;
    x= randn(6,1);
    w = 0.2+rand(6,1);
    p=randperm(6);
    param.g_d=p;
    param.g_t=[2 2 2];
    x1=[x(p(1)),x(p(2))];
    x2=[x(p(3)),x(p(4))];
    x3=[x(p(5)),x(p(6))];
    
    w1=[w(p(1)),w(p(2))];
    w2=[w(p(3)),w(p(4))];
    w3=[w(p(5)),w(p(6))];
    
    
    [ Kw1, ny1] = find_mg( x1',2*lambda,w1' );
    tg1 = 2*lambda/(1+2*lambda*Kw1)*ny1;
    
    [ Kw2, ny2] = find_mg( x2',2*lambda,w2' );
    tg2 = 2*lambda/(1+2*lambda*Kw2)*ny2;
    
    [ Kw3, ny3] = find_mg( x3',2*lambda,w3' );
    tg3 = 2*lambda/(1+2*lambda*Kw3)*ny3;
    
    n1=softC(x1,tg1);
    n2=softC(x2,tg2);
    n3=softC(x3,tg3);
    n=[n1,n2,n3];
    
    nx=n(back_perm(p'))';
    
    param.verbose=0;
    param.weights = w;
    
    %prox_l21(x,lambda,param);
    if (norm(nx- prox_l12(x,lambda,param))<epsilon)
        fprintf('Test L12 prox - Group same size - ok\n')
        b=0;
    else
        fprintf('Test L12 prox - Group same size - Error !!!!!!!!!!!!!!\n')
        b=1;
    end
     
    x= rand(6,1);
    p=randperm(6);
    w = 0.2+rand(6,1);
    
    param.g_d=p;
    param.g_t=[2 4];
    x1=[x(p(1)),x(p(2))];
    x2=[x(p(3)),x(p(4)),x(p(5)),x(p(6))];
    
    w1=[w(p(1)),w(p(2))];
    w2=[w(p(3)),w(p(4)),w(p(5)),w(p(6))];
    
    [ Kw1, ny1] = find_mg( x1',2*lambda,w1' );
    tg1 = 2*lambda/(1+2*lambda*Kw1)*ny1;
    
    [ Kw2, ny2] = find_mg( x2',2*lambda,w2' );
    tg2 = 2*lambda/(1+2*lambda*Kw2)*ny2;
    
    
    n1=softC(x1,tg1);
    n2=softC(x2,tg2);
    
    n=[n1,n2];
    nx=n(back_perm(p'))';
    
    param.weights = w;

    if (norm(nx- prox_l12(x,lambda,param))<epsilon)
        fprintf('Test L12 prox - Group diff size - ok\n')
        b=0+b;
    else
        fprintf('Test L12 prox - Group diff size - Error !!!!!!!!!!!!!!\n')
        b=1+b;
    end
 
end


function [b] = test_l12norm_matrix(epsilon)

    x= rand(6,10);
    
    nx=norm(sum(abs(x),2),2);
    
    

    if abs(nx-norm_l12(x))<epsilon
        fprintf('Test L12 norm - matrix - ok\n')
        b=0;
    else
        fprintf('Test L12 norm - matrix - Error !!!!!!!!!!!!!!\n')
        b=1;
    end


end


function [b] = test_l12prox_matrix(epsilon)

    b = 0;

    softC = @(z, T) sign(z).*max(abs(z)-T, 0);   
    lambda=0.5;

  
    param.verbose = 0;
    
    N =2;
    
    K = 10;
    
    x = rand(N,K);
    
    n = zeros(N,K);
    for ii = 1:N
        [ Kw, ny] = find_mg( x(ii,:)',2*lambda);
        tg1 = 2*lambda/(1+2*lambda*Kw)*ny;
        n(ii,:) = softC(x(ii,:)',tg1)';
    end
    

    
    n2 = prox_l12(x,lambda,param);



    if norm(n(:) - n2(:)) <epsilon
        fprintf('Test L12 prox - matrix 1 - ok\n')
    else
        fprintf('Test L12 prox - matrix 1 - Error !!!!!!!!!!!!!!\n')
        b=b+1;
    end
    
    
    
    x = 5*rand(100,1);


    param.verbose = 0;

  
    
    n2 = prox_l12(x,lambda,param);
    n = prox_l2(x,lambda,param);


    if norm(n(:) - n2(:)) <epsilon
        fprintf('Test L12 prox - matrix 2 - ok\n')
    else
        fprintf('Test L12 prox - matrix 2 - Error !!!!!!!!!!!!!!\n')
        b=b+1;
    end

end



function [b] = test_find_mg(epsilon)

    b = 0;

    x = randn(10);
    w = ones(10);
    [ ~, ~, ~, n,~] = find_mg( x,1,w );

    n2 = sort(abs(x),'descend');
    
    if norm(abs(n(:)) - n2(:)) <epsilon
        fprintf('Test find mg 1 - ok\n')
    else
        fprintf('Test find mg 1 - Error !!!!!!!!!!!!!!\n')
        b=b+1;
    end

    N = 100;
    x = 2 * randn(N,1);
    w = 0.2 + rand(N,1);
    
    lambda = 1.2;
    
    [ Kw, ny, Mg, zo, wo] = find_mg( x,lambda,w );
    
    ro = abs(zo./wo);
    
    temp1 = lambda * sum( wo(1:(Mg+1)).^2 .* ( ro(1:(Mg+1)) - ro(Mg+1) ) );
    temp2 = lambda * sum( wo(1:Mg).^2 .* ( ro(1:(Mg)) - ro(Mg) ) );

    
    if temp1 - ro(Mg+1)>=0
        fprintf('Test find mg 2 - ok\n')
    else
        fprintf('Test find mg 2 - Error !!!!!!!!!!!!!!\n')
        temp1 - ro(Mg)
        b=b+1;
    end
    
    
    
    if temp2 - ro(Mg)< 0
        fprintf('Test find mg 3 - ok\n')
    else
        fprintf('Test find mg 3 - Error !!!!!!!!!!!!!!\n')
        b=b+1;
    end
    
    ny2 = norm(zo(1:Mg),1);
    
    if norm(ny(:) - ny2(:)) <epsilon
        fprintf('Test find mg 4 - ok\n')
    else
        fprintf('Test find mg 4 - Error !!!!!!!!!!!!!!\n')
        b=b+1;
    end
    
    Kw2 = sum(wo(1:Mg).^2);

    if norm(Kw(:) - Kw2(:)) <epsilon
        fprintf('Test find mg 5 - ok\n')
    else
        fprintf('Test find mg 5 - Error !!!!!!!!!!!!!!\n')
        b=b+1;
    end
    
    
    N = 100;
    K = 10;
    
    x = 2 * randn(N,K);
    w = 0.2 + rand(N,K);
    
    lambda = 1.2;
    
    [ Kw, ny, Mg] = find_mg( x,lambda,w );
    
    Kw2 = zeros(K,1);
    ny2 = zeros(K,1);
    Mg2 = zeros(K,1);
    for ii = 1:K 
        [ Kw2(ii), ny2(ii), Mg2(ii)] = find_mg( x(:,ii),lambda,w(:,ii) );
    end
    
    if norm(Kw(:) - Kw2(:)) <epsilon
        fprintf('Test find mg 6 - ok\n')
    else
        fprintf('Test find mg 6 - Error !!!!!!!!!!!!!!\n')
        b=b+1;
    end
    
    if norm(ny(:) - ny2(:)) <epsilon
        fprintf('Test find mg 7 - ok\n')
    else
        fprintf('Test find mg 7 - Error !!!!!!!!!!!!!!\n')
        b=b+1;
    end
    
     
    if norm(Mg(:) - Mg2(:)) <epsilon
        fprintf('Test find mg 8 - ok\n')
    else
        fprintf('Test find mg 8 - Error !!!!!!!!!!!!!!\n')
        b=b+1;
    end
end
