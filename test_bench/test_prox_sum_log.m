function errors = test_prox_sum_log()

errors = 0;


errors = errors + errors_positive_input() + errors_negative_input();

end


function errors = errors_positive_input()

x = rand(100);

gamma = 0.1;
params.verbose = 2;

y = prox_sum_log(x, gamma, params);

errors = 0;
if not( all( y > x ))
    errors = errors + 1;
    fprintf('(WRONG)    prox(x) > x\n');
end
if not( all( y >= sqrt(4*gamma)/2))
    errors = errors + 1;
    fprintf('(WRONG)    prox(x) > sqrt(4*gamma)/2\n');
end


if errors == 0
    fprintf('(OK) elementwise magnitudes seem correct\n');
end

end

function errors = errors_negative_input()

x = randn(100);

gamma = 0.1;
params.verbose = 2;

y = prox_sum_log(x, gamma, params);

errors = 0;
if not( all( y > x ))
    errors = errors + 1;
    fprintf('(WRONG)    prox(x) > x\n');
end


if errors == 0
    fprintf('(OK) elementwise magnitudes seem correct\n');
end

end
