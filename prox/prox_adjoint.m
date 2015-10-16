function [sol]  = prox_adjoint(x, gamma, f)
%PROX_ADJOINT Proximal operator of the adjoint function of f
%   Usage:   sol=prox_adjoint(x, gamma, f);
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         f     : Function
%   Output parameters:
%         sol   : Solution.
%         infos : Structure summarizing informations at convergence
%
%   ` prox_adjoint( x,gamma,f )` solves:
%
%   .. sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * f*
%
%   .. math::  sol = \min_{z} \frac{1}{2} \|x - z\|_2^2 + \gamma  f^*
%
%   where $f^*$ is the adjoint of $f$. This problem is solved thanks to the
%   Moreau's identity.
%
%   Warning: f needs to be a proper convex lower semi continuous function.
%


% Author: Nathanael Perraudin
% Date: 31 May 2013
% Testing: test_prox_functions

if gamma == 0
    warning(['gamma = 0. This is problably not correct.' ...
        ' We replace it by eps to keep going.']);
    gamma = eps;
end

% Optional input arguments
if nargin < 3,
    error('Two few input arguments!');
end

% if nargout>1
%     [sol_a, info] = gamma*f.prox(x/gamma,1/gamma);
%     info.adjoint = 1;
% else
sol_a = gamma * f.prox(x / gamma, 1 / gamma);
% end
sol = x - sol_a;

end

