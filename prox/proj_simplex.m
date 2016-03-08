function [sol, info] = proj_simplex(x, ~, param)
%PROJ_SIMPLEX projection onto the simplex sum(x) = c, x >= 0, c = scalar
%   Usage:  sol = proj_simplex(x, ~, param)
%           [sol, info] = proj_simplex(x, ~, param)
%
%   Input parameters:
%         x     : Input signal.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         info : Structure summarizing informations at convergence
%
%   `proj_simplex(x,~,param)` solves:
%
%   .. sol = argmin_{z} ||x - z||_2^2   s.t.  sum(x) = c, x >= 0
%
%   .. math::  sol = \min_z ||x - z||_2^2 \hspace{1cm} s.t. \hspace{1cm} \sum_i x_i = c, x_i \geq 0
%
%   If x is a matrix, the projection is done for each column separately,
%   or across param.dim if specified.
%
%   *param* is a Matlab structure containing the following fields:
%
%   * *param.c* : scalar (default: 1).
%
%   * *param.dim* : dimension of summation if x is a matrix (default: 1)
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence (default: 0)
%
%   
%   info is a Matlab structure containing the following fields:
%
%   * *info.time* : Time of execution of the function in sec.
%
%   * *info.final_eval* : Final evaluation of the function
%
%
%
%
%   Rem: The input "~" is useless but needed for compatibility reasons.
%
%   See also:  proj_linear_eq, proj_b1
%
%   References: duchi2008efficient

% TODO: implement "figure 2" algorithm of same paper!! (difficult to
% vectorize!)
%
% Author: Vassilis Kalofolias
% Date: February 2016
% Testing: test_proj_simplex

% Start the time counter
t1 = tic;

if nargin < 3
    param = struct;
end


if ~isfield(param, 'c'), param.c = 1; end
if ~isfield(param, 'dim')
    if isvector(x)
        % 1 for column, 2 for row
        param.dim = 1 + isrow(x);
    else
        param.dim = 1;
    end
end
if ~isfield(param, 'verbose'), param.verbose = 0; end


[m, n] = size(x);

% implementing the algorithm of figure 1 of [1]
x_s = sort(x, param.dim, 'descend');
if param.dim == 1
    one_over_j = 1./(1:m)';
else
    one_over_j = 1./(1:n);
end
    
Thetas = bsxfun(@times, one_over_j, (cumsum(x_s, param.dim)-param.c)  );
% find number of points that are below 0
rho = sum((x_s - Thetas) > 0, param.dim);

if param.dim == 1
    theta = Thetas(sub2ind([m, n], rho, 1:n));
else
    theta = Thetas(sub2ind([m, n], (1:m)', rho));
end
sol = max(0, bsxfun(@plus, x, -theta));

% info about algorithm
info.algo = mfilename;
info.iter = 1;
info.final_eval = 0;
info.crit = '';
info.time = toc(t1);

if param.verbose >= 1
    fprintf('  Proj. simplex:  %e  seconds\n', info.time);
end

end


