function [ errors ] = test_proxlp( )
%TEST_PROX_LP This function test the function prox_lp
errors=0;


errors=errors+test_prox_lp_p2();

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -- test prox lp for p=2--  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [errors]=test_prox_lp_p2()
N=100;
x= randn(N,1);
y= randn(N,1);
param.y=y;
param.verbose=0;
param.p=2;

A=rand(N);
param.tight=0;
param.nu=norm(A)^2;
param.A=@(x) A*x;
param.At=@(x) A'*x;
param.tol = eps(1000);
param.maxit=1000;
%param.pcg = 0;
sol=prox_l2(x,1,param);
sol2=prox_lp(x,1,param);

    if norm(sol2(:)-sol(:))/norm(sol(:))<1e-2
       fprintf('Test prox Lp - simple 2 - ok\n')
       errors=0;
    else
       fprintf('Test prox Lp - simple 2 - Error !!!!!!!!!!!!!!\n')
       errors=1;
       norm(sol2(:)-sol(:))/norm(sol(:))
    end

end

