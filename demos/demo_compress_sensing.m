%DEMO_COMPRESS_SENSING Compress sensing example using forward backward algorithm
%
%   We present a compress sensing example solved with the forward backward
%   solver.
%   The problem can be expressed as this
%
%   ..   argmin ||Ax-b||^2 + tau ||x||_1
%
%   .. math:: \operatorname{arg\,min}_x \|Ax-b\|^2 + \tau \|x\|_{1}
%  
%   Where b are the measurements and A the measurement matrix.
%
%   We set 
%
%   * $f_1(x)=||x||_1$
%     We define the prox of $f_1$ as: 
%
%     .. prox_{f1,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma ||z||_1
%
%     .. math:: prox_{f1,\gamma} (z) = \operatorname{arg\,min}_{x} \frac{1}{2} \|x-z\|_2^2  +  \gamma \|z\|_{1}
%   
%     This function is simply a soft thresholding.
%
%   * $f_2(x)=||Ax-b||_2^2$
%     We define the gradient as: 
%
%     .. grad_f(x) = 2 * A^*(Ax-b)
%
%     .. math:: \nabla_f(x) = 2 A^*(x-b)
%
%     A is the measurement matrix (random Gaussian distribution)
%
%   The number of measurements $M$ is computed with respect of the size of the
%   signal $N$ and the sparsity level $K$:
%   
%   .. M=K*max(4,ceil(log(N)))   
%
%   .. math:: M=K \max\left(4,\text{ceil}(\log(N))\right)  
%
%   With this number of measurements, the algorithm is supposed to perform
%   very often always a perfect reconstruction. This plot is automatically
%   generated; let's hope it will be the case.
%
%   Results
%   -------
%
%   .. figure::
%
%      Results of the algorithm
%
%      This figure shows the original signal and the reconstruction done 
%      thanks to the algorithm and the measurements. The number of
%      measurements is M=900, the length of the signal N=5000 and K=100. This is
%      equivalent to a compression ratio of 5.55.
%
%   References: combettes2007douglas combettes2011proximal

% Author: Nathanael Perraudin, Gilles Puy
% Date: Nov 2012 


%% Initialisation

clear all;
close all;

% Loading toolbox
init_unlocbox();

verbose = 2; % verbosity level

%% Creation of the problem

tau = 1;                    % regularization parameter for the problem

N = 5000;                   % Size of the signal
K = 100;                    % Sparsity level
R = max(4, ceil(log(N)));   % Constant 
fprintf('The compression ratio is: %g\n',N/(R*K));


% Mesurements matrix
A = randn(R * K, N);

% Create a K sparse signal
x = zeros(N, 1); 
I = randperm(N);
x(I(1:K)) = randn(K, 1);
x = x / norm(x);

% Measurements
y = A * x;

%% Defining proximal operators

% setting the function f2 
f2.grad = @(x) 2*A'*(A*x-y);
f2.eval = @(x) norm(A*x-y)^2;
f2.beta = 2 * norm(A)^2;


% setting the function f1
param_l1.verbose = verbose -1;
param_l1.tight = 1;

f1.prox=@(x, T) prox_l1(x, T*tau, param_l1);
f1.eval=@(x) tau*norm(x,1);   


%% solving the problem

% setting different parameter for the simulation
param_solver.verbose = verbose; % display parameter
param_solver.maxit = 300;       % maximum iteration
param_solver.tol = 1e-4;        % tolerance to stop iterating
param_solver.method = 'FISTA';  % desired method for solving the problem

% solving the problem
sol = solvep(zeros(N,1), {f1, f2}, param_solver);

%% displaying the result
% figure;
plot(1:N, x, 'o', 1:N, sol, 'xr');
legend('Original signal', 'Reconstructed signal');
    
%% Closing the toolbox
close_unlocbox();
