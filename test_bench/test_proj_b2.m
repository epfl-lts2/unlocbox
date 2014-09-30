function [ errors ] = test_proj_b2( )
%TEST_PROX_L2 This function test the function prox_l2
errors=0;
tol=10e-5;

errors=errors+test_proj_b2_simple(tol);
errors=errors+test_proj_b2_FISTA(tol);


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -- test proj b2 simple --  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [errors]=test_proj_b2_simple(tol)

nx=30;
ny=15;
x= randn(nx,1);
y= randn(ny,1);
A=rand(ny,nx);
param.A=@(x) A*x;
param.At=@(x) A'*x;
param.nu=norm(A)^2;
param.y=y;
param.verbose=0;
param.epsilon=5;
param.tight=0;
param.maxit=1000;
param.tol=tol/10;

sol1=proj_b2(x,0,param);

param.method='ISTA';
sol2=proj_b2(x,0,param);
    if norm(sol1-sol2)/norm(sol1)<tol
       fprintf('Test proj b2 - simple - ok\n')
       errors=0;
    else
       fprintf('Test proj b2 - simple - Error !!!!!!!!!!!!!!\n')
       errors=1;
       norm(sol1-sol2)/norm(sol1)
    end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -- test proj b2 FISTA   --  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [errors]=test_proj_b2_FISTA(tol)

nx=30;
ny=15;
x= randn(nx,1);
y= randn(ny,1);
A=rand(ny,nx);
param.A=@(x) A*x;
param.At=@(x) A'*x;
param.nu=norm(A)^2;
param.y=y;
param.verbose=0;
param.epsilon=0.1;
param.tight=0;
param.maxit=1000;

sol1=fast_proj_b2(x,0,param);

param.method='FISTA';
sol2=proj_b2(x,0,param);
    if norm(sol1-sol2)/norm(sol1)<tol
       fprintf('Test proj b2 - FISTA - ok\n')
       errors=0;
    else
       fprintf('Test proj b2 - FISTA - Error !!!!!!!!!!!!!!\n')
       errors=1;
       norm(sol1-sol2)/norm(sol1)
    end

end

