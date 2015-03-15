%DEMO_RLR  Example of use of the RLR solver
%
%   We present an example of the RLR solver through an image reconstruction problem.
%   The problem can be expressed as this
%
%   ..   argmin ||Ax-b||^2 + tau*||x||_TV
%
%   .. math:: arg \min_x \|Ax-b\|^2 + \tau \|x\|_{TV}
%  
%   Where b is the degraded image, I the identity and A an operator representing the mask.
%
%   We set 
%
%   * $f_1(x)=||x||_{TV}$
%     We define the prox of $f_1$ as: 
%
%     .. prox_{f1,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma ||z||_TV
%
%     .. math:: prox_{f1,\gamma} (z) = arg \min_{x} \frac{1}{2} \|x-z\|_2^2  +  \gamma \|z\|_{TV}
%
%   Results
%   -------
%
%   .. figure::
%
%      Original image
%
%      This figure shows the original Lena image. 
%
%   .. figure::
%
%      Depleted image
%
%      This figure shows the image after the application of the mask. Note
%      that 50% of the pixels have been removed.
%
%   .. figure::
%
%      Reconstruted image
%
%      This figure shows the reconstructed image thanks to the algorithm.
%   
%
%   References: combettes2011proximal

% Author: Nathanael Perraudin
% Date: November 2012
%

%% Initialisation

clear all;
close all;

% Loading toolbox
init_unlocbox();

verbose = 2;    % verbosity level

%% Defining the problem

tau = 0.1; % Regularization parameter for the problem

% Original image
im_original = lena();

% Creating the problem
A=rand(size(im_original));
A=(A>0.5);
% Depleted image
b=A.*im_original;

% Defining the adjoint operator
A_op = @(x) A.*x;
At_op = @(x) A.*x;

%% Define the proximal operator

% setting the function f1 (norm TV)
param_tv.verbose = verbose - 1;
param_tv.maxit = 50;

f.prox=@(x, T) prox_tv(x, T*tau, param_tv);
f.eval=@(x) tau*norm_tv(x);   

%% Solving the problem

% setting different parameters for the simulation
param.verbose = verbose;    % display parameter
param.maxit = 40;           % maximum iteration
param.epsilon = 10e-5;      % tolerance to stop iterating
param.gamma = 0.5;          % stepsize (beta is equal to 2)
param.method = 'FISTA';     % desired method for solving the problem

% solving the problem
sol=rlr(b,f,A_op,At_op,param);

%% displaying the result
imagesc_gray(im_original, 1, 'Original image');
imagesc_gray(b, 2, 'Depleted image');
imagesc_gray(sol, 3, 'Reconstructed image');
    
%% Closing the toolbox
close_unlocbox();
