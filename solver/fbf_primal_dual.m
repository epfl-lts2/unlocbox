function [sol, info] = fbf_primal_dual(x_0,f1, f2, f3, param)
%FBF_PRIMAL_DUAL forward backward forward primal dual
%   Usage: sol = fbf_primal_dual(x_0,f1, f2, f3, param);
%          sol = fbf_primal_dual(x_0,f1, f2, f3);
%          [sol,info,objective] = fbf_primal_dual(...);
%
%   Input parameters:
%         x_0   : Starting point of the algorithm
%         f1    : First function to minimize
%         f2    : Second function to minimize
%         f3    : Third function to minimize
%         param : Optional parameters
%   Output parameters:
%         sol   : Solution
%         info  : Structure summarizing informations at convergence
%
%   `fbf_primal_dual` (using forward backward forward based primal dual)
%   solves:
%
%   .. sol = argmin f1(x) + f2(Lx) + f3(x)
%
%   .. math::  sol = \min_x f_1(x) + f_2( L x) + f_3(x)
%   
%   where  $x$ is the optimization variable with $f_1$ or $f_3$ a smooth
%   function and $L$ a linear operator. $f_1$ and $f_3$ are defined like
%   other traditional functions.
%
%   Note that *f2* is a structure of a functions with:
%
%   * `f2.eval(x_i)` : an operator to evaluate the function
%   * `f2.prox(x_i, gamma)` : an operator to evaluate the prox of the function
%
%   Optionally you can define
%
%   * `f2.L`  : linear operator, matrix or operator (default identity)
%   * `f2.Lt` : adjoint of linear operator, matrix or operator (default identity)
%   * `f2.norm_L` : bound on the norm of the operator L (default: 1), i.e.
%
%     .. ` ||L x||^2 <= nu * ||x||^2 
%
%     .. math::  \|L x\|^2 \leq \nu \|x\|^2 
%
%
%   *param* a Matlab structure containing solver paremeters. See the
%   function |solvep| for more information. Additionally it contains those
%   aditional fields:  
%
%   * *param.tol* : is stopping criterion for the loop. The algorithm stops if
%
%     ..  max_i ||  y_i(t) - y_i(t-1) ||  / ||y(t) ||< tol,
%      
%     .. math:: \max_i \frac{ \| y_i(t) - y_i(t-1)\| }{ \|y_i(t)\|} < tol,
%
%     where  $y_i(t)$ are the dual variable of function *i* at itertion *t*
%     by default, `tol=10e-4`.
%
%       Warning! This stopping criterion is different from other solvers!
%
%   * *param.mu* : parameter mu of paper [1]
%   * *param.epsilon*:   parameter epsilon of paper [1]
%   * *param.normalized_timestep*: from 0 to 1, mapping to [epsilon,
%     (1-epsilon)/mu]
%
%   See also: solvep fb_based_primal_dual
%
%   Demos:  demo_fbb_primal_dual
%
%   References: komodakis2014playing
 
% Author: Vassilis Kalofolias, Nathanael Perraudin
% Date: June 2015
% Testing: test_solvers

param.algo = 'FBF_PRIMAL_DUAL';
[sol, info] = solvep(x_0,{f1, f2, f3},param);

end