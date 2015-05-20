function [sol, info,objective] = fb_based_primal_dual(x_0,f1, f2, f3, param)
%FB_BASED_PRIMAL_DUAL forward backward based primal dual
%   Usage: sol = fb_based_primal_dual(x_0,f1,f2, f3,param);
%          sol = fb_based_primal_dual(x_0,f1,f2,f3);
%          [sol,info,objective] = fb_based_primal_dual(...);
%
%   Input parameters:
%         x_0   : Starting point of the algorithm
%         f1    : First function to minimize
%         f2    : Second function to minimize
%         f3    : Third function to minimize
%         param : Optional parameter
%   Output parameters:
%         sol   : Solution
%         info  : Structure summarizing informations at convergence
%         objective: vector (evaluation of the objectiv function each iteration)
%
%   `admm` (using alternating-direction method of multipliers) solves:
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
%
%   *param* a Matlab structure containing solver paremeters. See the
%   function |solvep| for more information. Additionally it contains those
%   aditional fields:  
%
%   * *param.tol* : is stop criterion for the loop. The algorithm stops if
%
%     ..  max_i ||  y_i(t) - y_i(t-1) ||  / ||y(t) ||< tol,
%      
%     .. math:: \max_i \frac{ \| y_i(t) - y_i(t-1)\| }{ \|y_i(t)\|} < tol,
%
%     where  $y_i(t)$ are the dual variable of function *i* at itertion *t*
%     by default, `tol=10e-4`.
%
%       Warning! This stopping criterion is different from other solver!
%
%   * *param.sigma* : second timestep. By default equal to *param.gamma*.
%
%   * *param.rescale* : Use the rescale version of the algorithm (default 0)
%
%   See also: solvep sdmm admm
%
%   Demos:  demo_admm
%
%   References: komodakis2014playing
 
% Author: Nathanael Perraudin
% Date: 2 May 2015
% Testing: test_solvers

param.algo = 'FB_BASED_PRIMAL_DUAL';
[sol, info,objective] = solvep(x_0,{f1, f2, f3},param);

end