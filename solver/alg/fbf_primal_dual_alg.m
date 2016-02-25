function s = fbf_primal_dual_alg()
%FBF_PRIMAL_DUAL Forward-Backward-Forward primal dual algorithm
%   Usage : param.algo = fbf_primal_dual();
%
%   This function returns a structure containing the algorithm. You can
%   lauch your personal algorithm with the following::
%
%           param.algo = fbf_primal_dual_alg();
%           sol = solvep(x0, {f1, f2, f3}, param);
%

%   code authors: Vassilis Kalofolias, Nathanael Perraudin
%   date: June 2015
% This function returns a structure with 4 fields:
% 1) The name of the solver. This is used to select the solvers.
s.name = 'FBF_PRIMAL_DUAL';
% 2) A method to initialize the solver (called at the beginning)
s.initialize = @(x_0, fg, Fp, param) ...
    fbf_primal_dual_initialize(x_0,fg,Fp,param);
% 3) The algorithm itself (called at each iterations)
s.algorithm = @(x_0, fg, Fp, sol, s, param) ...
    fbf_primal_dual_algorithm(fg, Fp, sol, s, param);
% 4) Post process method (called at the end)
s.finalize = @(x_0, fg, Fp, sol, s, param) sol;
% The variables here are
%   x_0 : The starting point
%   fg  : A single smooth function
%         (if fg.beta == 0, no smooth function is specified)
%   Fp  : The non smooth functions (a cell array of structure)
%   param: The structure of optional parameters
%       .mu:        parameter mu of paper [1]
%       .epsilon:   parameter epsilon of paper [1]
%       .normalized_timestep: from 0 to 1, mapping to [epsilon,
%                               (1-epsilon)/mu]
%   s   : Intern variables or the algorithm
%   sol : Current solution
%
%
% [1] playing with duality, primal dual ... Komodakis, Pesquet.
%
% see also: fbf_primal_dual
end




function [sol, s, param] = fbf_primal_dual_initialize(x_0, fg, Fp, param)

if (numel(Fp)>2)
    error('This solver needs at maximum 2 non-smooth functions')
end

% add a dummy non-smooth function if none given
if (numel(Fp)==1)
    Fp{2}.prox = @(x) x;
    Fp{2}.eval = eps;
end

% struct keeping the information of the different parts of the
% objective function
%
% ind: indexing so that the second one (the last one?) has the linear
%       operator L
% L:    linear operator used to go to the dual space
% Lt:   adjoint of L operator to go from dual to primal space
s = struct;

% Reorder functions so that second one is the one with transformation L
if isfield(Fp{1},'L')
    s.ind = [2,1];
    L = Fp{1}.L;
    Lt = Fp{1}.Lt;
elseif isfield(Fp{2},'L')
    s.ind = [1,2];
    L = Fp{2}.L;
    Lt = Fp{2}.Lt;
    % add dummy if no L used
else
    L = @(x) x;
    Lt = @(x) x;
    s.ind = [1,2];
end


if isfield(Fp{s.ind(2)}, 'norm_L')
    s.norm_L = Fp{s.ind(2)}.norm_L;
else
    s.norm_L = 1;
    if isL
        warning('You should give f.norm_L = ||L||^2. Setting it to 1!');
    end
end


% compute a timestep
beta = fg.beta;

if not(isfield(param, 'mu'))
    s.mu = beta + s.norm_L;   % TODO: check that norm is not squared
else
    s.mu = param.mu;
end
if not(isfield(param, 'epsilon'))
    s.epsilon = 0;       % in (0, 1/(1+mu) )
else
    s.epsilon = param.epsilon;
end

% timestep
if not(isfield(param, 'normalized_timestep'))
    param.normalized_timestep = 0.5;
end

s.tau =  lin_map(param.normalized_timestep, [s.epsilon, (1-s.epsilon)/s.mu], [0, 1]);    % in [epsilon, (1-epsilon)/mu]


% All internal variables are stored into the structure s
s.OpL = L;
s.OpLt = Lt;
sol = x_0;                  % primal variable

s.dual_var = L(sol);        % dual variable
s.g_prox_adjoint = @(x,T) prox_adjoint(x,T,Fp{s.ind(2)});


end

function [sol, s] = fbf_primal_dual_algorithm(fg, Fp, sol, s, param)

if (numel(Fp)==1)
    Fp{2}.prox = @(x) x;
    Fp{2}.eval = eps;
end

%TODO: set s.tau here

Y_n = sol - s.tau * (fg.grad(sol) + s.OpLt(s.dual_var));
y_n = s.dual_var + s.tau * (s.OpL(sol));
P_n = Fp{s.ind(1)}.prox(Y_n, s.tau);
p_n = s.g_prox_adjoint(y_n, s.tau);
Q_n = P_n - s.tau * (fg.grad(P_n) + s.OpLt(p_n));
q_n = p_n + s.tau * (s.OpL(P_n));

sol = sol - Y_n + Q_n;
s.dual_var = s.dual_var - y_n + q_n;
end