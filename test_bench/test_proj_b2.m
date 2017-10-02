function [ errors ] = test_proj_b2( )
%TEST_PROX_L2 This function test the function prox_l2
errors=0;
tol=10e-5;

errors=errors+test_proj_b2_simple(tol);
errors=errors+test_proj_b2_FISTA(tol);

gsp_reset_seed(1)
errors=errors+test1();



end


function errors = test1()

errors = 0;

y = randn(100,1);
z = zeros(100,1);
epsilon = 1.3;

paramp.y = z;
paramp.epsilon = epsilon;
paramp.tight = 1;

x = proj_b2(y,1,paramp);
errors = errors + assert_test(norm(x-z),epsilon,1e-8,'Test proj_b2 1');

z = randn(100,1);
paramp.y = z;

x = proj_b2(y,1,paramp);
errors = errors + assert_test(norm(x-z),epsilon,1e-8,'Test proj_b2 2');

A = dctmtx(100);
paramp.A = @(x) A*x;
paramp.At = @(x) A'*x;

x = proj_b2(y,1,paramp);
errors = errors + assert_test(norm(A*x-z),epsilon,1e-8,'Test proj_b2 3');


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

