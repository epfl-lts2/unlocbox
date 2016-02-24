function [ errors ] = test_prox_l2( )
%TEST_PROX_L2 This function test the function prox_l2
errors=0;
tol=10e-10;

errors=errors+test_prox_l2_simple(tol);
errors=errors+test_prox_l2_weights(tol);
errors=errors+test_prox_l2_y(tol);
errors=errors+test_prox_l2_iter1(tol);
errors=errors+test_prox_l2_iter2(tol);
errors=errors+test_prox_l2_iter3(tol);

errors=errors+test_prox_l2_tightT(tol);


errors=errors+test_prox_l2_op1(tol);
errors=errors+test_prox_l2_op2(tol);
errors=errors+test_prox_l2_op3(tol);
errors=errors+test_prox_l2_op4(tol);
errors=errors+test_prox_l2_op5(tol);

errors=errors+test_prox_l2_pcg(tol);
errors=errors+test_prox_l2_comp_pcg(tol);


errors=errors+test_prox_l2_multidim(tol);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -- test prox l2 simple --  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [errors]=test_prox_l2_simple(tol)
x= randn(100,1);
sol=x/3;
param.verbose=0;
    if norm(prox_l2(x,1,param)-sol)<tol
       fprintf('Test prox L2 - simple - ok\n')
       errors=0;
    else
       fprintf('Test prox L2 - simple - Error !!!!!!!!!!!!!!\n')
       errors=1;
    end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -- test prox l2 y       --  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [errors]=test_prox_l2_y(tol)
x= randn(100,1);
w= rand(100,1)+0.1;
y= randn(100,1);
gamma=0.5;
sol=(x+2*y*gamma.*w.^2)./(gamma*2*w.^2+1);
param.verbose=0;
param.weights=w;
param.y=y;
param.tol=eps(10);
    if norm(prox_l2(x,gamma,param)-sol)<1e-6
       fprintf('Test prox L2 - y - ok\n')
       errors=0;
    else
       fprintf('Test prox L2 - y - Error !!!!!!!!!!!!!!\n')
       param.verbose=2;
       norm(prox_l2(x,gamma,param)-sol)
       errors=1;
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -- test prox l2 y       --  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [errors]=test_prox_l2_tightT(tol)
x= randn(100,1);
y= randn(100,1);
w= rand(100,1);

M = rand(100,1);
M = double(M>0.5);
gamma=0.5;

param.A = @(x) M.*x;
param.At= @(x) M.*x;
param.verbose=0;
param.weights=w;
param.y=y;
param.tight = 0;
param.tol=eps(10);

sol = prox_l2(x,gamma,param);

param.tightT = 1;
    if norm(prox_l2(x,gamma,param)-sol)<1e-6
       fprintf('Test prox L2 - TightT - ok\n')
       errors=0;
    else
       fprintf('Test prox L2 - TightT - Error !!!!!!!!!!!!!!\n')
       param.verbose=2;
       norm(prox_l2(x,gamma,param)-sol)
       errors=1;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -- test prox l2 weights --  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [errors]=test_prox_l2_weights(tol)
x= randn(100,1);
w= rand(100,1)+0.1;

sol=x./(0.5*2*w.^2+1);
param.tight=0;
param.verbose=0;
param.weights=w;
param.tol=eps(10);
    if norm(prox_l2(x,0.5,param)-sol)<1e-6
       fprintf('Test prox L2 - weights - ok\n')
       errors=0;
    else
       fprintf('Test prox L2 - weights - Error !!!!!!!!!!!!!!\n')
       param.verbose=2;
       norm(prox_l2(x,0.5,param)-sol)
       errors=1;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -- test prox l2 iter 1 --  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [errors]=test_prox_l2_iter1(tol)
x= randn(100,1);
gamma=0.5;
sol=x/(1+2*gamma);
param.verbose=0;
param.tight=0;
param.tol=eps(10);
    if norm(prox_l2(x,gamma,param)-sol)<tol
       fprintf('Test prox L2 - iter 1 - ok\n')
       errors=0;
      
    else
       fprintf('Test prox L2 - iter 1 - Error !!!!!!!!!!!!!!\n')
       errors=1;
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -- test prox l2 iter 2  --  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [errors]=test_prox_l2_iter2(tol)
x= randn(100,1);
w= rand(100,1);
gamma=1;
sol=x./(gamma*2*w.^2+1);
param.verbose=0;
param.weights=w;
param.tight=0;
param.tol=eps(10);
    if norm(prox_l2(x,gamma,param)-sol)<tol
       fprintf('Test prox L2 - iter 2 - ok\n')
       errors=0;
    else
       fprintf('Test prox L2 - iter 2 - Error !!!!!!!!!!!!!!\n')
       errors=1;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -- test prox l2 iter 3  --  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [errors]=test_prox_l2_iter3(tol)
x= randn(100,1);
w= rand(100,1);
y= randn(100,1);
gamma=0.5;
sol=(x+2*y*gamma.*w.^2)./(gamma*2*w.^2+1);
param.verbose=0;
param.weights=w;
param.y=y;
param.tight=0;
param.tol=eps(10);
    if norm(prox_l2(x,gamma,param)-sol)<tol
       fprintf('Test prox L2 - iter 3 - ok\n')
       errors=0;
    else
       fprintf('Test prox L2 - iter 3 - Error !!!!!!!!!!!!!!\n')
       errors=1;
    end
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -- test prox l2 op 1    --  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [errors]=test_prox_l2_op1(tol)
L=100;
x= randn(L,1);
w= rand(L,1);
y= randn(L,1);
gamma=0.5;
sol=(gamma*2*ifft(diag(w.^2)*fft(eye(L)))+eye(L))^(-1)*(x+2*gamma.*sqrt(L).*ifft(y.*w.^2));
param.verbose=0;
param.weights=w;
param.y=y;
param.tol=eps(10);
param.A = @(x) 1/sqrt(100)*fft(x);
param.At = @(x) sqrt(100)*ifft(x);
    if norm(prox_l2(x,gamma,param)-sol)<tol
       fprintf('Test prox L2 - op 1 - ok\n')
       errors=0;
    else
       fprintf('Test prox L2 - op 1 - Error !!!!!!!!!!!!!!\n')
       errors=1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -- test prox l2 op 2    --  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [errors]=test_prox_l2_op2(tol)
L=100;
x= randn(L,1);
w= rand(L,1);
y= randn(L,1);
gamma=0.5;
sol=(gamma*2*ifft(diag(w.^2)*fft(eye(L)))+eye(L))^(-1)*(x+2*gamma.*sqrt(100).*ifft(y.*w.^2));
param.verbose=0;
param.weights=w;
param.y=y;
param.tight=0;
param.tol=eps(10);
param.A = @(x) 1/sqrt(100)*fft(x);
param.At = @(x) sqrt(100)*ifft(x);
    if norm(prox_l2(x,gamma,param)-sol)<tol
       fprintf('Test prox L2 - op 2 - ok\n')
       errors=0;
    else
       fprintf('Test prox L2 - op 2 - Error !!!!!!!!!!!!!!\n')
       errors=1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -- test prox l2 op 3    --  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [errors]=test_prox_l2_op3(tol)
L=100;
x= randn(L,1);

gamma=0.5;

% create a tight frame
a=25;
M=50;
g=gabwin({'gauss',1},a,M,L);
gt=gabtight(g,a,M,L);
gt=gt/norm(gt);
Fr=frame('dgt',gt,a,M);
F=frsynmatrix(Fr,L);
[A,B]=gabframebounds(gt,a,M,L);


w= ones(M*L/a,1);
% set the parameters
param.A = @(x) F'*x;
param.At = @(x) F*x;
param.nu=A;
y= randn(M/a*L,1);


param.verbose=0;
param.weights=w;
param.y=y;
param.tol=0;
param.tight=1;
sol=(2*gamma*F*diag(w.^2)*F'+eye(L))^(-1)*(x+2*gamma.*F*(w.^2.*y));

    if norm(prox_l2(x,gamma,param)-sol)<tol
       fprintf('Test prox L2 - op 3 - ok\n')
       errors=0;
    else
       fprintf('Test prox L2 - op  3 - Error !!!!!!!!!!!!!!\n')
       errors=1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -- test prox l2 op 4    --  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [errors]=test_prox_l2_op4(tol)
L=100;
x= randn(L,1);

gamma=0.5;

% create a tight frame
a=25;
M=50;
g=gabwin({'gauss',1},a,M,L);
gt=gabtight(g,a,M,L);
gt=gt/norm(gt);
Fr=frame('dgt',gt,a,M);
F=frsynmatrix(Fr,L);
[A,B]=gabframebounds(gt,a,M,L);


w= rand(M*L/a,1);
% set the parameters
param.A = @(x) F'*x;
param.At = @(x) F*x;
param.nu=A;
y= randn(M/a*L,1);


param.verbose=0;
param.weights=w;
param.y=y;
param.tol=eps(10);
param.tight=0;
sol=(2*gamma*F*diag(w.^2)*F'+eye(L))^(-1)*(x+2*gamma.*F*(w.^2.*y));

    if norm(prox_l2(x,gamma,param)-sol)<tol
       fprintf('Test prox L2 - op 4 - ok\n')
       errors=0;
    else
       fprintf('Test prox L2 - op  4 - Error !!!!!!!!!!!!!!\n')
       errors=1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -- test prox l2 op 5    --  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [errors]=test_prox_l2_op5(tol)
L=64;
x= randn(L,1);
w= 1+rand(L,1);
%w=ones(L,1);
% set the parameters
y= randn(L,1);
%y=zeros(L,1);

param.verbose=0;
param.weights=w;
param.y=y;
param.tol = eps(100);



f.prox=@(x,gamma) prox_l2(x,gamma,param);
f.eval=@(x) norm(x-y)^2;
f2.prox=@(x,T) x;
f2.eval=@(x) eps;

% param_p.abs_tol=1;
param_p.tol=0.5*tol;
param_p.verbose=0;
sol=ppxa(x,{f,f2},param_p);

    if abs(f.eval(sol))<tol
       fprintf('Test prox L2 - op 5 - ok\n')
       errors=0;
    else
       fprintf('Test prox L2 - op  5 - Error !!!!!!!!!!!!!!\n')
       errors=1;
       abs(f.eval(sol))
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -- test prox l2 pcg   --  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [errors]=test_prox_l2_pcg(tol)
L=100;
x= randn(L,1);

gamma=0.5;

% create a tight frame
a=25;
M=50;
g=gabwin({'gauss',1},a,M,L);
gt=gabtight(g,a,M,L);
gt=gt/norm(gt);
Fr=frame('dgt',gt,a,M);
F=frsynmatrix(Fr,L);
[A,B]=gabframebounds(gt,a,M,L);


w= ones(M*L/a,1);
% w=ones(size(x));

% F= eye(L);

% set the parameters
param.A = @(x) F'*x;
param.At = @(x) F*x;
param.nu=B;
y= randn(M/a*L,1);
% y = zeros(size(x));

param.verbose=0;
param.weights=w;
param.y=y;
param.tol=0.1*tol;
param.tight=0;
param.pcg=1;
sol=(2*gamma*F*diag(w.^2)*F'+eye(L))^(-1)*(x+2*gamma.*F*(w.^2.*y));

    if norm(prox_l2(x,gamma,param)-sol)<tol
       fprintf('Test prox L2 - pcg - ok\n')
       errors=0;
    else
       fprintf('Test prox L2 - pcg - Error !!!!!!!!!!!!!!\n')
       errors=1;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -- test comp pcg   --  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [errors]=test_prox_l2_comp_pcg(tol)
L=100;
N=200;
x= randn(L,1);

gamma=0.5;

F=randn(N,L);

A=norm(F)^2;


w= rand(N,1);
% set the parameters
param.A = @(x) F*x;
param.At = @(x) F'*x;
param.nu=A;
y= randn(N,1);


param.verbose=1;
param.weights=w;
param.y=y;
param.tol=100*eps;
param.tight=0;
param.pcg=1;
param.maxit=10000;
[sol1,infos1] = prox_l2(x,gamma,param);

param.pcg=0;
[sol2,infos2] = prox_l2(x,gamma,param);


    if norm(sol1-sol2)<1e-4
       fprintf('Test prox L2 - comp pcg - ok\n')
       fprintf('   Time with gradient descent: %g\n',infos2.time);
       fprintf('   Time with pcg             : %g\n',infos1.time);
       errors=0;
    else
       fprintf('Test prox L2 - comp pcg - Error !!!!!!!!!!!!!!\n')
       norm(sol1-sol2)
       errors=1;
    end
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -- test prox L2 multidim --  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [errors]=test_prox_l2_multidim(tol)
L=100;
k = 5;
x= randn(L,k);

gamma=0.5;

% create a tight frame
a=25;
M=50;
g=gabwin({'gauss',1},a,M,L);
gt=gabtight(g,a,M,L);
gt=gt/norm(gt);
Fr=frame('dgt',gt,a,M);
F=frsynmatrix(Fr,L);
[A,B]=gabframebounds(gt,a,M,L);


w= ones(M*L/a,k);
% w=ones(size(x));

% F= eye(L);

% set the parameters
param.A = @(x) F'*x;
param.At = @(x) F*x;
param.nu=B;
y= randn(M/a*L,k);
% y = zeros(size(x));

param.verbose=0;
param.weights=w;
param.y=y;
param.tol=0.1*tol;
param.tight=0;
param.pcg=1;
sol1=prox_l2(x,gamma,param);

sol2 = zeros(L,k);
for ii = 1:k
    param.y = y(:,ii);
    param.weights = w(:,ii);
    sol2(:,ii)=prox_l2(x(:,ii),gamma,param);
end


    if norm(sol1(:)-sol2(:))/norm(sol1(:))<tol
       fprintf('Test prox L2 - pcg - ok\n')
       errors=0;
    else
       fprintf('Test prox L2 - pcg - Error !!!!!!!!!!!!!!\n')
       errors=1;
       norm(sol1(:)-sol2(:))/norm(sol1(:))
    end
end


