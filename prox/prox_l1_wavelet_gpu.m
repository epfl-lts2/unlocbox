function [sol,infos] = prox_l1_wavelet_gpu(x, gamma, param)
%PROX_L1_WAVELET_GPU Proximal operator with L1 norm
%   Usage:  sol=prox_l1(x, gamma)
%           sol=prox_l1(x, gamma, param)
%           [sol, infos]=prox_l1(x, gamma, param)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         infos : Structure summarizing informations at convergence
%
%   This function is experimental!
%
%   `prox_l1(x, gamma, param)` solves:
%
%   .. sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * ||A z||_1
%
%   .. math::  sol = \min_{z} \frac{1}{2} \|x - z\|_2^2 + \gamma  \|A z\|_1
%
%   param is a Matlab structure containing the following fields:
%
%   * *param.A* : Forward operator (default: Id).
%
%   * *param.At* : Adjoint operator (default: Id).
%
%   * *param.tight* : 1 if A is a tight frame or 0 if not (default = 1)
%
%   * *param.nu* : bound on the norm of the operator A (default: 1), i.e.
%
%     .. ` ||A x||^2 <= nu * ||x||^2 
%
%     .. math::  \|A x\|^2 \leq \nu  \|x\|^2 
%   
%   * *param.tol* : is stop criterion for the loop. The algorithm stops if
%
%     ..  (  n(t) - n(t-1) )  / n(t) < tol,
%      
%     .. math:: \frac{  n(t) - n(t-1) }{ n(t)} < tol,
%
%     where  $n(t) = f(x)+ 0.5 \|x-z\|_2^2$ is the objective function at iteration *t*
%     by default, `tol=10e-4`.
%
%   * *param.maxit* : max. nb. of iterations (default: 200).
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%   * *param.weights* : weights for a weighted L1-norm (default = 1)
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
%
%   References: fadili2009monotone beck2009fast


% Author:  Eyal Hirsch, Nathanael Perraudin
% Date: Jan 2013
%

% Start the time counter
t1 = tic;

% Optional input arguments
if nargin<3, param=struct; end

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'tight'), param.tight = 1; end
if ~isfield(param, 'nu'), param.nu = 1; end
if ~isfield(param, 'tol'), param.tol = 1e-3; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'At'), param.At = @(x) x; end
if ~isfield(param, 'A'), param.A = @(x) x; end
if ~isfield(param, 'weights'), param.weights = 1; end
if ~isfield(param, 'pos'), param.pos = 0; end

% test the parameters
gamma=test_gamma(gamma);
param.weights=test_weights(param.weights);

% Useful functions

    % This soft thresholding function only support real signal
    % soft = @(z, T) sign(z).*max(abs(z)-T, 0);

    % This soft thresholding function support complex number
    soft = @(z, T) max(abs(z)-T,0)./(max(abs(z)-T,0)+T).*z;


% Projection
if param.tight && ~param.pos % TIGHT FRAME CASE
    
    temp = param.A(x);
    sol = x + 1/param.nu * param.At(soft(temp, ...
        gamma*param.nu*param.weights)-temp);
    crit = 'REL_OBJ'; iter = 0;
    dummy = param.A(sol);
    norm_l1 = sum(param.weights(:).*abs(dummy(:)));

else
    % NON TIGHT FRAME CASE OR CONSTRAINT INVOLVED
    % Initializations - The original code takes significant time and
    % should either be optimized (as below) or moved to the GPU.
    % Original code:
    %   u_l1 = zeros(size(param.Psit(x)));
    %   sol = x - param.Psi(u_l1);
    % Optimized code:
    u_l1 = zeros(size(x));
    sol = x;

    
    
    % Soft-thresholding
    % Init
    if param.verbose > 1
        fprintf('  Proximal l1 operator:\n');
    end
    
    if sum(imag(sol(:)))
        warning('Data converted to real values\n');       
    end
    sol = real(sol);
    x = real(x);

    Depth = 1;
    Width  = size(sol, 1);
    Height = size(sol, 2);
    qmfSize = size(param.h, 2);
    result=single(ones(size(sol)));
    result_ptr = libpointer('singlePtr', result);
    x_ptr = libpointer('singlePtr', single(x));
    sol_ptr = libpointer('singlePtr', single(sol));
    u_l1_ptr = libpointer('singlePtr', single(u_l1));
    qmf_ptr = libpointer('singlePtr', single(param.h));
    gpuResults=zeros(4, 1);
    gpuResults_ptr = libpointer('singlePtr', single(gpuResults));
    
    binaryName=get_gpu_binary_name();
    result = calllib(binaryName,'ProxL1', result_ptr, x_ptr, sol_ptr, u_l1_ptr, gpuResults_ptr, ...
               Width, Height, Depth, gamma, param.weights, param.nu, param.tol, param.maxit, param.pos, param.L, qmf_ptr, qmfSize);
    sol = reshape(result_ptr.Value, [Width Height Depth]);

    % Parse the return values from the GPU implementation.
    gpuResults = reshape(gpuResults_ptr.Value, [size(gpuResults, 1) size(gpuResults, 2)]);
    norm_l1=gpuResults(1);
    rel_l1=gpuResults(2);
    iter=gpuResults(3);
    if (rel_l1 < param.tol)
        crit = 'REL_OB';
    elseif iter >= param.maxit
        crit = 'MAX_IT';
    end
end

% Log after the projection onto the L2-ball
if param.verbose >= 1
    fprintf(['  GPU prox_L1: ||A x||_1 = %e,', ...
        ' %s, iter = %i\n'], norm_l1, crit, iter);
end



infos.algo=mfilename;
infos.iter=iter;
infos.final_eval=norm_l1;
infos.crit=crit;
infos.time=toc(t1);

end
