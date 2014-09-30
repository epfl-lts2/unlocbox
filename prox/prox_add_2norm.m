function [ sol]  = prox_add_2norm( x,gamma,param )
%PROX_ADD_2NORM Proximal operator with an additional quadratic term
%   Usage:   sol=prox_add_2norm(x, gamma, param);
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         f     : Function
%   Output parameters:
%         sol   : Solution.
%         infos : Structure summarizing informations at convergence
%
%   ` prox_add_2norm( x,gamma,param )` solves:
%
%   .. sol = argmin_{z} 0.5*||x - z||_2^2 + 0.5*||y - z||_2^2 + gamma * f
%
%   .. math::  sol = \min_{z} \frac{1}{2} \|x - z\|_2^2 +\frac{1}{2} \|y - z\|_2^2+ \gamma  f
%
%   This problem can be solved because we have the nice relationship
%
%   .. 0.5*||x - z||_2^2 + 0.5*||y - z||_2^2 = || (x+y)/2 - z||_2^2 
%   ..                                         + 0.25 *||y - x||_2^2
%
%   .. math:: \frac{1}{2} \|x - z\|_2^2 +\frac{1}{2} \|y - z\|_2^2 =  \|\frac{x+y}{2} - z\|_2^2 +\frac{1}{4} \|y - x\|_2^2
%
%   This function can be used to reduce the number of functional to
%   minimize and thus to employ faster solvers.
%
%   *param* is a Matlab structure containing the following fields:
%   
%   * *param.y* : a vector of the same size as x
%
%   * *param.f* : a structure containing the function f
%

% Author: Nathanael Perraudin
% Date: 1 February 2014
%



% Optional input arguments
if nargin<3,     
    error('Two few input arguments!');
end

if ~isfield(param, 'y')
   error('Please specify param.y')
end
if ~isfield(param, 'f')
   error('Please specify param.f')
end

sol = param.f.prox((x+param.y)/2,0.5*gamma);

end

