function [ errors ] = test_lp( )
errors=0;

gsp_reset_seed(1);
errors=errors+test_eq();
errors=errors+test_ineq();




end

function [errors]=test_eq()
    
    errors = 0;

    A = rand(10,50);
    x = rand(50,1);
    y = rand(10,1);
    
    param.pinvA = pinv(A);
    param.A = A;
    param.y = y;
    param.method = 'exact';
    
    sol = proj_linear_eq(x,0,param);
    
    if norm(A*sol-y)/norm(sol)<1e-10
        fprintf('  Test proj_linear_eq  OK\n')
    else
        fprintf('  Test proj_linear_eq Pas OK!!!!!!!!!!!!!!!!\n')
        norm(A*sol-y)/norm(sol)<1e-10
        errors= errors +1;
    end
    
    param.method = 'proj_b2';
    param.At = A';
    param.maxit = 2000;
    param.verbose = 0;
    param.tol = 1e-10;
    param.nu = norm(A)^2;
    param.tight = 0;
    sol2 = proj_linear_eq(x,0,param);
    if norm(sol-sol2)/norm(sol)<1e-5
        fprintf('  Test proj_linear_eq proj_b2 OK\n')
    else
        fprintf('  Test proj_linear_eq proj_b2 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(sol-sol2)/norm(sol)
        errors= errors +1;
    end
    
    param.method = 'primal_dual';
    param.At = A';
    param.maxit = 2000;
    param.verbose = 0;
    param.tol = 1e-10;
    param.nu = norm(A)^2;
    param.tight = 0;
    sol2 = proj_linear_eq(x,0,param);
    if norm(sol-sol2)/norm(sol)<1e-5
        fprintf('  Test proj_linear_eq primal_dual OK\n')
    else
        fprintf('  Test proj_linear_eq primal_dual Pas OK!!!!!!!!!!!!!!!!\n')
        norm(sol-sol2)/norm(sol)
        errors= errors +1;
    end
end

function [errors]=test_ineq()
    
    errors = 0;

    A = randn(10,50);
    x = randn(50,1);
    y = -randn(10,1);
    
    param.pinvA = pinv(A);
    param.A = A;
    param.y = y;
    param.method = 'quadprog';
    
    sol = proj_linear_ineq(x,0,param);
    if ~sum(A*sol-y>1e-5)
        fprintf('  Test proj_linear_ineq  OK\n')
    else
        fprintf('  Test proj_linear_ineq Pas OK!!!!!!!!!!!!!!!!\n')
        A*sol-y
        errors= errors +1;
    end
    
    param.method = 'iterative';
    param.At = A';
    param.maxit = 2000;
    param.verbose = 1;
    param.tol = 1e-10;
    param.nu = norm(A)^2;
    param.tight = 0;
    sol2 = proj_linear_ineq(x,0,param);
    
    if norm(sol-sol2)/norm(sol)<1e-5
        fprintf('  Test proj_linear_ineq iterative OK\n')
    else
        fprintf('  Test proj_linear_ineq iterative Pas OK!!!!!!!!!!!!!!!!\n')
        norm(sol-sol2)/norm(sol)
        errors= errors +1;
    end
end
