%DEMO_BPDN Deconvolution demonstration (Debluring)
%
%   Here we try to deblur an image through a deconvolution problem. The
%   convolution operator is the blur
%   The problem can be expressed as this
%
%   ..   argmin  ||Ax-b||^2 + tau*||H(x)||_1
%
%   .. math:: arg \min_x \|Ax-b\|^2 + \tau \|H(x)\|_{1}
%  
%   Where b is the degraded image, I the identity and A an operator representing the blur.
%
%   H is a linear operator projecting the signal in a sparse
%   representation. Here we worked with wavelet. 
%
%   Warning! Note that this demo require the LTFAT to work.
%
%
%   Results
%   -------
%
%   .. figure::
%
%      Original image
%
%      This figure shows the original lena image. 
%
%   .. figure::
%
%      Depleted image
%
%      This figure shows the image after the application of the blur.
%
%   .. figure::
%
%      Reconstructed image
%
%      This figure shows the reconstructed image thanks to the algorithm.
%

% Author: Nathanael Perraudin
% Date: 14 May 2013
%

%% Initialisation
clear;
close all;

% Loading toolbox
init_unlocbox();

verbose = 1;    % Verbosity level

%% Creation of the problem

% Original image
im_original = cameraman();

% Creation of a Gaussian round blur
sigma = 0.1; % width of the blur
[x, y] = meshgrid(linspace(-1, 1, length(im_original)));
r = x.^2 + y.^2;
G = exp( -r / (2 * sigma^2) );

% Blur operators
A = @(x) real(ifft2(fftshift(G).* (fft2(x))));
At = @(x) real(ifft2(conj(fftshift(G)).*(fft2(x))));
    
% Depleted image
b = A(im_original);

%% Defining proximal operators

% Prior assumption about the image (sparse in Wavelet)
L=8;
Psi = @(x) fwt2(x,'db1',L);
Psit = @(x) ifwt2(x,'db1',L);
%%

% setting different parameter for the simulation
% General parameter
param.verbose = verbose;    % display parameter
param.maxit  =50;     % maximum iteration
param.tol = 1e-5;     % tolerance to stop iterating
param.gamma = 0.03;   % Timestep (very important parameter)

% Paramter for the meaurements
epsilon = 0;          % fidelity (radius of the B2 ball)
param.nu_b2 = 1;      % Lipshitz constant of A
param.tight_b2 = 0;   % If A is tight
param.maxit_b2 = 100; % Maximum number of iteration for the projection 
                      % onto the ball b2.
param.tol_b2 = 1e-5;  % Tolerance to stop iterating


% Parameter for the L1 prior assumption
param.maxit_l1 = 50;  % Maximum number of itertaion for the minimization 
                      % of the L1 norm.
param.tight_l1 = 1;   % If Psi is tight
param.nu_l1 = 1;      % See documentation of solve_bpdn 

      

% solving the problem
sol = solve_bpdn(b, epsilon, A, At, Psi, Psit, param);
    

%% displaying the result
imagesc_gray(im_original, 1, 'Original image');
imagesc_gray(b, 2, 'Depleted image');
imagesc_gray(sol, 3, 'Reconstructed image');
    

%% Closing the toolbox
close_unlocbox();

