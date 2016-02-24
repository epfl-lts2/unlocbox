% test for proj_simplex.m
%
% code author: Vassilis Kalofolias
% date: Feb 2016

function num_errors = test_proj_simplex()

num_errors = 0;

tol = 1e-10;

num_errors = num_errors + test_vector(tol) + test_matrix(tol);



end

function num_errors = test_vector(tol)
    
    num_errors = 0;

    
    %% row vector
    x = randn(100, 1);
    y = proj_simplex(x);
    
    if any(y < -tol)
        num_errors = num_errors + 1;
    end
    
    if abs((sum(y) - 1)) > tol
        num_errors = num_errors + 1;
    end
    
    y = proj_simplex(x');
    
    if any(y < -tol)
        num_errors = num_errors + 1;
    end
    
    if abs((sum(y) - 1)) > tol
        num_errors = num_errors + 1;
    end
    
    if num_errors == 0
        fprintf('  Test proj_simplex for vectors OK\n')
    else
        fprintf('  Test proj_simplex for vectors NOT OK!!!!!!!!\n')
    end

    
end


function num_errors = test_matrix(tol)
    
    num_errors = 0;

    
    X = randn(100, 10);
    
    params.c = 3;       % default is 1
    Y1 = proj_simplex(X, [], params);
    
    % if I transpose and do the same across second dimension and transpose
    % again I should get the same result!
    params.dim = 2; % default is 1
    Y2 = proj_simplex(X', [], params)';
    
    
      
    if any(Y1(:) < -tol)
        num_errors = num_errors + 1;
    end
      
    if norm(Y1 - Y2, 'fro') / norm(Y1, 'fro') > tol
        num_errors = num_errors + 1;
    end
    
    if norm((sum(Y1) - params.c)) > tol
        num_errors = num_errors + 1;
    end
    
    if num_errors == 0
        fprintf('  Test proj_simplex for matrices OK\n')
    else
        fprintf('  Test proj_simplex for matrices NOT OK!!!!!!!!\n')
    end

    
end



