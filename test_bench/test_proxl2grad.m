function [ errors ] = test_proxl2grad( )
%TEST_PROXL2GRAD test the function prox_l2grad
errors=0;
tol=10e-10;

errors=errors+test_proxl2grad_1dim(tol);
errors=errors+test_proxl2grad_1dimop(tol);
errors=errors+test_proxl2grad_1dimopreverse(tol);
errors=errors+test_proxl2grad_1dimopmulticol(tol);

errors=errors+test_proxl2grad_2dimop(tol);




end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -- test prox proxl2grad 1dim simple--  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [errors]=test_proxl2grad_1dim(tol)
N=100;
x= randn(N,1);

sol=prox_l2grad_old(x,1);
sol2=prox_l2grad(x,1);

    if norm(sol2-sol)<tol
       fprintf('Test prox l2 grad - 1 dim - ok\n')
       errors=0;
    else
       fprintf('Test prox l2grad - 1 dim - Error !!!!!!!!!!!!!!\n')
       errors=1;
       norm(sol2-sol)
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -- test prox proxl2grad 1dim operator--  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [errors]=test_proxl2grad_1dimop(tol)
N=100;
x= randn(N,1);
param.A=@(x) 1/sqrt(N)*fft(x);
param.At=@(x) sqrt(N)*ifft(x);

sol2=prox_l2grad(x,1,param);
sol=prox_l2gradfourier(x,1);

    if norm(sol2-sol)<tol
       fprintf('Test prox l2 grad - 1 dim op - ok\n')
       errors=0;
    else
       fprintf('Test prox l2grad - 1 dim op - Error !!!!!!!!!!!!!!\n')
       errors=1;
       norm(sol2-sol)
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -- test prox proxl2grad 1dim reverse--  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [errors]=test_proxl2grad_1dimopreverse(tol)
N=100;
x= randn(N,1);
param.A=@(x) 1/sqrt(N)*fft(x);
param.At=@(x) sqrt(N)*ifft(x);

sol2=prox_l2grad(x,1,param);
sol=prox_l2grad(x',1,param);

    if norm(sol2-sol')<tol
       fprintf('Test prox l2 grad - 1 dim op reverse - ok\n')
       errors=0;
    else
       fprintf('Test prox l2grad - 1 dim op reverse - Error !!!!!!!!!!!!!!\n')
       errors=1;
       norm(sol2-sol)
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -- test prox proxl2grad 1dim operator--  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [errors]=test_proxl2grad_1dimopmulticol(tol)
N=100;
x= randn(N,30);
param.A=@(x) 1/sqrt(N)*fft(x);
param.At=@(x) sqrt(N)*ifft(x);

sol2=prox_l2grad(x,1,param);
sol=prox_l2gradfourier(x,1);

    if norm(sol2-sol)<tol
       fprintf('Test prox l2 grad - 1 dim op mulitcol - ok\n')
       errors=0;
    else
       fprintf('Test prox l2grad - 1 dim op multicol - Error !!!!!!!!!!!!!!\n')
       errors=1;
       norm(sol2-sol)
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -- test prox proxl2grad 2dim operator--  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [errors]=test_proxl2grad_2dimop(tol)
N=10;
L=20;
x= randn(N,L);
param.A=@(x) 1/sqrt(N)*1/sqrt(L)*fft2(x);
param.At=@(x) sqrt(L)*sqrt(N)*ifft2(x);
param.d2 = 1;

sol2=prox_l2grad(x,1,param);
sol=prox_l2gradfourier(x,1,param);

    if norm(sol2-sol)<tol
       fprintf('Test prox l2 grad - 2 dim op - ok\n')
       errors=0;
    else
       fprintf('Test prox l2grad - 2 dim op - Error !!!!!!!!!!!!!!\n')
       errors=1;
       norm(sol2-sol)
    end

end

