function [ errors ] = test_proj_b1( )
%TEST_PROJ_B1 test prox_b1 and related function
errors=0;

errors=errors+test_simple(eps(100));




end

function [errors]=test_simple(tol)
    
    errors = 0;

    x = randn(100,1);
    tau = norm(x,1)/1.5;
    param.verbose = 0;
    param.epsilon = tau;
    
    p3 = oneProjector(x,1,tau);
    p2 = proj_b1(x,0,param);
    
    
    if norm(p3-p2)/norm(p3)<tol
        fprintf('  Test proj_b1 simple 1 OK\n')
    else
        fprintf('  Test proj_b1 simple Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
    x = randn(100,1);
    weight = randn(100,1);
    tau = norm(weight.* x,1)/1.5;
    param.verbose = 0;
    param.epsilon = tau;
    param.weight = weight;
    
    p3 = oneProjector(x,weight,tau);
    p2 = proj_b1(x,0,param);
    
    
    if norm(p3-p2)/norm(p3)<tol
        fprintf('  Test proj_b1 simple 2 OK\n')
    else
        fprintf('  Test proj_b1 simple 2 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    n3 = norm(weight.*p2,1);
    if  n3<tau*1.01;
        fprintf('  Test proj_b1 simple 3 OK\n')
    else
        fprintf('  Test proj_b1 simple 3 Pas OK!!!!!!!!!!!!!!!!\n')
        errors= errors +1;
    end
    
    
    x = randn(100,1);
    weight = randn(100,1);
    tau = norm(weight.* x,1);
    param.verbose = 0;
    param.epsilon = tau*1.1;
    param.weight = weight;
    
    p3 = oneProjector(x,weight,tau);
    p2 = proj_b1(x,0,param);
    
    if norm(p3-p2)/norm(p3)<tol
        fprintf('  Test proj_b1 simple 4 OK\n')
    else
        fprintf('  Test proj_b1 simple 4 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
    if norm(p2-x)/norm(x)<tol
        fprintf('  Test proj_b1 simple 5 OK\n')
    else
        fprintf('  Test proj_b1 simple 5 Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p2-x)/norm(x)
        errors= errors +1;
    end
end


