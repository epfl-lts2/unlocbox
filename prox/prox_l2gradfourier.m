function [sol,info] = prox_l2gradfourier(x, gamma, param)
%PROX_L2gradfourier Proximal operator of the 2 norm of the gradient in the Fourier domain
%   Usage:  sol=prox_l2gradfourier(x, gamma)
%           sol=prox_l2gradfourier(x, gamma, param)
%           [sol, info]=prox_l2gradfourier(x, gamma, param)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         info  : Structure summarizing informations at convergence
%
%   This function compute the 1 dimensional proximal operator of x. For
%   matrices, the function is applied to each column. The parameter
%   param.d2 allows the user to use the 2 dimentional gradient.
%
%   Warning: the signal should not be centered. Indice 1 for abscissa 0.
%
%   `prox_l2gradfourier(x, gamma, param)` solves:
%
%   .. sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * ||grad(Fz)||_2^2
%
%   .. math::  sol = \min_{z} \frac{1}{2} \|x - z\|_2^2 + \gamma  \|\nabla \mathcal{F}z\|_2^2
%
%   param is a Matlab structure containing the following fields:
%
%   * *param.weights* : weights if you use a an array.
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%   * *param.deriveorder* : Order ot the derivative default 1
%
%   * *param.d2* : 2 dimentional gradient (default 0)
%
%
%   info is a Matlab structure containing the following fields:
%
%   * *info.algo* : Algorithm used
%
%   * *info.iter* : Number of iteration
%
%   * *info.time* : Time of exectution of the function in sec.
%
%   * *info.final_eval* : Final evaluation of the function
%
%   * *info.crit* : Stopping critterion used 
%
%
%   See also:  prox_l2 prox_l2grad prox_tv



% Author: Nathanael Perraudin
% Date: Jan 2013
% Testing: test_prox_l2grad

% Start the time counter
t1 = tic;

% Optional input arguments
if nargin<3, param=struct; end

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'weights'), param.weights = 1; end
if ~isfield(param, 'deriveorder'), param.deriveorder = 1; end
if ~isfield(param, 'd2'), param.d2 = 0; end



warning=0;
% test the parameters
if test_gamma(gamma)
    sol = x;
    info.algo=mfilename;
    info.iter=0;
    info.final_eval=0;
    info.crit='--';
    info.time=toc(t1);
    return; 
end

test_weights(param.weights,warning);


p=param.deriveorder;

% useful function
h=@(t) 1./(1+2*param.weights(:)*gamma*t).^p;

if param.d2
    % Make the frequencies vector
        % a) indices
        L=size(x,1);
        Q=size(x,2);
        l=(0:L-1)';
        q=(0:Q-1);
        % b) values
        eig_l = (2-2*cos(2*pi*l/L));
        eig_q = (2-2*cos(2*pi*q/Q));
        % c) Compute the radius for the kernel
        rho = repmat(eig_l,1,Q) + repmat(eig_q,L,1);
        sol=x.*repmat(h(rho),[1,1,size(x,3)]);
    
    [dx, dy] =  gradient_op(1/sqrt(L)*1/sqrt(Q)*fft2(sol));
    curr_norm = norm(dx.^2+dy.^2,'fro')^2;
elseif (size(x,1)==1) || (size(x,2)==1)
    L=size(x,1)*size(x,2);
    if size(x,1)>size(x,2)
        l=(0:L-1)';
    else
        l=(0:L-1);
    end
    lambda=2-2*cos(2*pi*l/L);        

    %filtering
    sol=x.*h(lambda);
    curr_norm = norm(gradient_op1d(1/sqrt(L)*fft(sol(:))),'fro')^2;
else
    L=size(x,1);
    l=(0:L-1)';
    lambda=2-2*cos(2*pi*l/L);
    sol = x.*repmat(h(lambda),1,size(x,2));
    curr_norm = norm(gradient_op1d(1/sqrt(L)*fft(sol)),'fro')^2;
end

% Summary
if param.verbose>=1
   fprintf('  Prox_l2gradfourier: 1 iteration, ||grad(Fx)||^2=%g\n',curr_norm);
end

% zero iteration
iter=1;
crit='--';
info.algo=mfilename;
info.iter=iter;
info.final_eval=curr_norm;
info.crit=crit;
info.time=toc(t1);

end



