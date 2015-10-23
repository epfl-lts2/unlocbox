 

% This function need to be rewritten

% scalar input and scalar constraints
x = -5
param.lower_lim = 0;
param.upper_lim = 2;
param

[sol, info] = proj_box(x, [], param)


% vector input and scalar constraints
x = -5:5
param.lower_lim = 0;
param.upper_lim = 2;
param

[sol, info] = proj_box(x, [], param)


% vector input and vector constraints
x = -5:5
param.lower_lim = [1 3 2 0 0 0 0 0 3 3 3];
param.upper_lim = 3*ones(1,length(x));
param

[sol, info] = proj_box(x, [], param)
