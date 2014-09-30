%DEMO_GRADIENT_DESCENT Simple 1 dimension deconvolution
%
%   In this demos, we use convex optimization to blindly deconvolve a 1
%   dimention signal. The method used is probably not a very efficient way
%   do perform deconvolution. However, under given conditions, it shoud
%   deconvolve perfectly the signal.
%
%   We try to solve:
%
%   ..   argmin  ||x*h-y||^2
%
%   .. math:: arg \min_x \|x*h-y\|^2 
%  
%   Where y is the convoluted signal and h the convolution kernel.
%
%   Let's write $x*h$ as Hx. We Set 
%
%   * $f(x)=||x*h-b||_2^2$
%     We define the gradient as: 
%
%     .. grad_f(x) = 2 H^*(Hx-y) 
%
%     .. math:: \nabla_f(x) = 2 H^*(Hx-y)
%
%   Results
%   -------
%
%   .. figure::
%
%      Convolution kernel
%
%      This figure shows the chosen konvolution kernel. It is a weighted
%      blur.
%
%   .. figure::
%
%      Original signal and convoluted signal
%
%      This figure shows modification of the signal after application of
%      the blur.
%
%   .. figure::
%
%      Original sisnal and solution to the gradient descent problem.
%
%      This figure shows the reconstructed signal thanks to the algorithm.


% Author: Nathanael Perraudin
% Date: sept 14 2012
%

%% Initialisation

clear all;
close all;

% Loading toolbox
init_unlocbox();

verbose = 2;    % verbosity level

%% Defining the problem

% Size of the problem
N=100;

% Creating signal.
s=rand(N,1);

% Creation of the blur
h=[1,3,4,1,1]';
%h=rand(N,1);

h=h/norm(h);

% Blur the signal
y=real(ifft(fft(h,N).*fft(s)));

%% Setting the gradient

% setting the function f
blurf=fft(h,N);
A=@(x) real(ifft(blurf.*(fft(x))));
At=@(x) real(ifft(conj(blurf).*(fft2(x))));

f.eval=@(x) norm(A(x)-y)^2;
f.grad=@(x) 2*At(A(x)-y);

%% Solving the problem

% setting different parameter for the simulation
param.verbose=0; % display parameter
param.maxit=1000; % maximum iteration
param.tol=10e-9; % tolerance to stop iterating
param.gamma=0.5/norm(diag(blurf))^2; % stepsize (beta is equal to 2)
param.method='FISTA'; % desired method for solving the problem

% solve the problem
sol=gradient_descent(y,f,param);

%% Display the results

% Displaying the blur
figure(1);
plot(1:N,real(ifft(fft(h,N))),'xb');
legend('Bluring kernel in time');


% Displaying depleted image
figure(2)
plot(1:N,s,'xb',1:N,y,'og');
legend('Original signal','Convoluted signal');

% displaying the result
figure(3)
plot(1:N,s,'xb',1:N,sol,'og');
legend('Original signal','Convluted signal');
fprintf('SNR of the reconstructed signal: %g dB\n',snr(s,sol));

