function [sol,infos] = prox_l2grad_old(x, gamma, param)
%PROX_L2grad_old Proximal operator of the 2 norm of the gradient in 1 dimension
%   Usage:  sol=prox_l2grad_old(x, gamma)
%           sol=prox_l2grad_old(x, gamma, param)
%           [sol, infos]=prox_l2grad_old(x, gamma, param)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         infos : Structure summarizing informations at convergence
%
%   This function compute the 1 dimensional proximal operator of x. For
%   matrices, the function is applied to each column. For N-D
%   arrays, the function operates on the first 
%   dimension.
%
%   `prox_l2grad(x, gamma, param)` solves:
%
%   .. sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * ||grad(z)||_2^2
%
%   .. math::  sol = \min_{z} \frac{1}{2} \|x - z\|_2^2 + \gamma  \|\nabla z\|_2^2
%
%   param is a Matlab structure containing the following fields:
%
%   * *param.abasis* : to use another basis than the DFT (default: 0). To be
%                     done -- Not working yet
%
%   * *param.weights* : weights if you use a an array.
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%   * *param.deriveorder* : Order ot the derivative default 1
%
%
%   infos is a Matlab structure containing the following fields:
%
%   * *infos.algo* : Algorithm used
%
%   * *param.iter* : Number of iteration
%
%   * *param.time* : Time of exectution of the function in sec.
%
%   * *param.final_eval* : Final evaluation of the function
%
%   * *param.crit* : Stopping critterion used 
%
%
%   See also:  proj_b1 prox_l1inf prox_l12 prox_tv


% Author: Nathanael Perraudin
% Date: Nov 2012
%

% Start the time counter
t1 = tic;

% Optional input arguments
if nargin<3, param=struct; end

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'weights'), param.weights = 1; end
if ~isfield(param, 'deriveorder'), param.deriveorder = 1; end
if ~isfield(param, 'abasis'), param.abasis = 0; end

warning=0;
% test the parameters
test_gamma(gamma);
test_weights(param.weights,warning);

if param.abasis
    error('This option is not currently supported, please contact the developer')
end

p=param.deriveorder;


% useful function
h=@(t) 1./(1+2*param.weights*gamma*t).^p;

% size of the signal




L=size(x,1);
l=(0:L-1)';
lambda=2-2*cos(2*pi*l/L);

y=fft(x);

sol=ifft(y.*h(lambda));


    

% one iteration
iter=1;

curr_norm=norm(gradient(sol),'fro')^2;

% Summary
if param.verbose>=1
   fprintf('  Prox_l2grad: %i iteration(s), ||grad(x)||^2=%g\n',iter,curr_norm);
end

crit='--';
iter=0;
infos.algo=mfilename;
infos.iter=iter;
infos.final_eval=curr_norm;
infos.crit=crit;
infos.time=toc(t1);

end



