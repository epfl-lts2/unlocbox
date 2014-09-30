function [ errors ] = test_op( )
%TEST_SOLVERS test differents solvers
errors=0;

errors=errors+test_laplacianx();
errors=errors+test_laplaciany();
errors=errors+test_laplacian();





end


function [errors]=test_laplacianx()
    errors = 0;
    
    N =100;
    
    x = rand(N,10);
    
    p2 = laplacianx_op(x);
    
    p3 = div_op1d(gradient_op1d(x));

    
    if norm(p3-p2)/norm(p3)<1e-6
        fprintf('  Test laplacianx 1 OK\n')
    else
        fprintf('  Test laplacianx Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end



function [errors]=test_laplaciany()
    errors = 0;
    
    N =100;
    
    x = rand(10,N);
    
    p2 = laplaciany_op(x);
    
    p3 = div_op1d(gradient_op1d(x'))';

    
    if norm(p3-p2)/norm(p3)<1e-6
        fprintf('  Test laplaciany 1 OK\n')
    else
        fprintf('  Test laplaciany Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3-p2)/norm(p3)
        errors= errors +1;
    end
    
end

function [errors]=test_laplacian()
    errors = 0;
    
    N =100;
    
    x = rand(N);
    
    p2 = laplacian_op(x);
    [a,b] = gradient_op(x);
    p3 = div_op(a,b);

    
    if norm(p3(:)-p2(:))/norm(p3(:))<1e-6
        fprintf('  Test laplacian OK\n')
    else
        fprintf('  Test laplacian Pas OK!!!!!!!!!!!!!!!!\n')
        norm(p3(:)-p2(:))/norm(p3(:))
        errors= errors +1;
    end
    
end
