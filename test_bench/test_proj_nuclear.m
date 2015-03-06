function errors = test_proj_nuclear()

errors = 0;


errors = errors + test_simple();
errors = errors + test_simple2();

end


function errors = test_simple()

errors = 0;
P = rand(10);

P = P+ P';

np = norm_nuclear(P);
param.verbose = 0;
param.epsilon = np*1.1;

res1 = proj_nuclearnorm(P,0,param);

if norm(res1 - P,'fro')<1e-3
    fprintf(' TEST_PROJ_NUCLEAR 1 OK\n')
else
    warning('Error in test_proj_nuclear 1!!!!!!!!!!')
    errors = errors + 1;
end


end


function errors = test_simple2()

errors = 0;
P = rand(100);

P = P+ P';

np = norm_nuclear(P);
param.verbose = 0;
param.epsilon = np/1.1;

res1 = proj_nuclearnorm(P,0,param);

if abs(norm_nuclear(res1)-param.epsilon)<1e-6
    fprintf(' TEST_PROJ_NUCLEAR 2 OK\n')
else
    warning('Error in test_proj_nuclear 2!!!!!!!!!!')
    errors = errors + 1;
end


end