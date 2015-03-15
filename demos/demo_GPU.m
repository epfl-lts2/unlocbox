%DEMO_GPU Demeonstration of the use of the GPU
%
%   Here we try to deblur an image through a deconvolution problem. The
%   convolution operator is the blur
%   The problem can be expressed as this
%
%   ..   argmin  ||Ax-b||^2 + tau*||x||_TV
%
%   .. math:: arg \min_x \|Ax-b\|^2 + \tau \|x)|_{TV}
%  
%   Where b is the degraded image and A an operator representing the blur.
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
%   * $f_2(x)=||Ax-b||_2^2$
%     We define the gradient as: 
%
%     .. grad_f(x) = 2 A^*(Ax-b)
%
%     .. math:: \nabla_f(x) = 2 A^*(Ax-b)
%
%   GPU
%   ---
%
%   We solve the problems two times. One with the GPU and and one without
%   and we compare also the time of execution. 
%
%   We use the GPU in two different way. First only the proximal tv
%   operator is computed with the GPU.
%
%   On my computer, I do not have a gain of time with the GPU. The
%   implementation might be disastrous. GPU acceleration on Matlab, we are
%   not there yet.
%
%   Results
%   -------
%
%   .. figure::
%
%      Original image
%
%      This figure shows the original128 lena image. 
%
%   .. figure::
%
%      Depleted image
%
%      This figure shows the image after the application of the blur.
%
%   .. figure::
%
%      Reconstructed image GPU
%
%      This figure shows the reconstructed image thanks to the algorithm.
%
%   .. figure::
%
%      Reconstructed image CPU
%
%      This figure shows the reconstructed image thanks to the algorithm.
%

% Author: Nathanael Perraudin, Gilles Puy
% Date: jan 23 2013
%

%% Initialisation

clear;
close all;

% Loading toolbox
init_unlocbox();

verbose = 1;    % verbosity level
global GLOBAL_useGPU;

%% Create a problem

% Original image
im_original = lena(); 

% Creating the problem
M=rand(size(im_original))<0.3;
A=@(x) M.*x;
At=A;

% Depleted image
b=A(im_original);

%% Setting proximal operator

tau = 0.1; % Regularization parameter for the problem

% setting the function f
param_tv.verbose = verbose - 1;
param_tv.maxit = 50;
param_tv.tol = 1e-5;

f.prox=@(x, T) prox_tv(x, T*tau, param_tv);
f.eval=@(x) tau * norm_tv(x);  

%% Solving the problem

% setting different parameters for the simulation
param_solver.verbose = verbose; % display parameter
param_solver.maxit = 20;        % maximum number of iterations
param_solver.tol = 10e-6;       % tolerance to stop iterating
param_solver.nu = 1;

%% solving the problem With the GPU
fprintf('Solve the problem with the GPU \n');
GLOBAL_useGPU = 1; %#ok<NASGU>
[solGPU,infos] = rlr(b,f,A,At,param_solver);
timeGPU = infos.time; 
fprintf('  Computation time with the GPU: %g\n',timeGPU);


%% solving the problem with CPU
fprintf('Solve the problem with the CPU only \n');
GLOBAL_useGPU = 0; 
[solCPU,infos] = rlr(b,f,A,At,param_solver);
timeCPU = infos.time; 
fprintf('  Computation time with the CPU: %g\n',timeCPU);

%% solving the problem all on the GPU
fprintf('Solve the problem all on the GPU only \n');
bgpu = gpuArray(b);
Mgpu = gpuArray(M);
Agpu=@(x) Mgpu.*x;
Atgpu=Agpu;
[solallGPU,infos] = rlr(bgpu,f,Agpu,Atgpu,param_solver);
timeallGPU = infos.time; 
solallGPU = gather(solallGPU);
fprintf('  Computation time all on the GPU: %g\n',timeallGPU);


%% Print result in the console
fprintf(' -------------- Summary -------------- \n\n');

fprintf('  Computation time with the GPU: %g\n',timeGPU);
fprintf('  Computation time with the CPU: %g\n',timeCPU);
fprintf('  Computation time all on the GPU: %g\n',timeallGPU);

fprintf('\n - Difference between the answer: \n');
fprintf('  SNR CPU / GPU : %d db\n', snr(solCPU,solGPU));
fprintf('  SNR CPU / all GPU : %d db\n', snr(solCPU,solallGPU));
fprintf('  SNR GPU / all GPU : %d db\n', snr(solGPU,solallGPU));

fprintf('\n ---------------- End ---------------- \n');

%% displaying the result
imagesc_gray(im_original, 1, 'Original image');
imagesc_gray(b, 2, 'Depleted image');
imagesc_gray(solGPU, 3, ...
    strcat('Reconstructed image with the GPU -- time:',num2str(timeGPU)));
imagesc_gray(solCPU, 4, ...
    strcat('Reconstructed image with the CPU -- time:',num2str(timeCPU)));
imagesc_gray(solallGPU, 5, ...
    strcat('Reconstructed image all on the GPU -- time:',num2str(timeallGPU)));

%% Closing the toolbox
close_unlocbox();
