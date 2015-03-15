%DEMO_COMPRESS_SENSING2 Compress sensing example using Douglas Rachford algorithm
%
%   We present a compress sensing example solved with the douglas rachford
%   solver.
%   The problem can be expressed as this
%
%   ..   argmin || x ||_1 s.t ||b-Ax||_2 < epsilon
%
%   .. math:: arg \min_x \| x\|_1 \hspace{1cm} such \hspace{0.25cm}  that \hspace{1cm} \|b-Ax\|_2 \leq \epsilon
%  
%   Where b are the measurements and A the measurement matrix.
%
%   We set 
%
%   * $f_1(x) = || x ||_1$
%     We define the prox of $f_1$ as: 
%
%     .. prox_{f1,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma ||z||_1
%
%     .. math:: prox_{f1,\gamma} (z) = arg \min_{x} \frac{1}{2} \|x-z\|_2^2  +  \gamma \|z\|_{1}
%   
%     This function is simply a soft thresholding.
%
%   * $f_2$ is the indicator function of the set S define by $||Ax-b||_2 < \epsilon$
%     We define the prox of $f_2$ as 
%
%     .. prox_{f2,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma i_S( x ),
%
%     .. math:: prox_{f2,\gamma} (z) = arg \min_{x} \frac{1}{2} \|x-z\|_2^2   + i_S(x) ,
%
%     with $i_S(x)$ is zero if x is in the set S and infinity otherwise.
%     This previous problem has an identical solution as:
%
%     .. argmin_{z} ||x - z||_2^2   s.t.  ||b - A z||_2 < epsilon
%
%     .. math:: arg \min_{z} \|x - z\|_2^2   \hspace{1cm} such \hspace{0.25cm} that \hspace{1cm} \|Az-b\|_2 \leq \epsilon
%
%     It is simply a projection on the B2-ball.
%     A is the measurement matrix (random Gaussian distribution)
%
%   The number of measurements $M$ is computed with respect of the size of the
%   signal $N$ and the sparsity level $K$:
%   
%   .. M=K*max(4,ceil(log(N)))   
%
%   .. math:: M=K \max\left( 4 , \text{ceil}(\log(N))\right)  
%
%   With this number of measurements, the algorithm is supposed to perform
%   very often always a perfect reconstruction. This plot is automatically
%   generated, let's hope it will be the case.
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
% Date: november 2012


%% Initialisation

clear;
close all;

% Loading toolbox
init_unlocbox();

verbose = 2; % verbosity level

%% Creation of the problem

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

% Define the prox of f2 see the function proj_B2 for more help
operatorA = @(x) A * x;
operatorAt = @(x) A' * x;
epsilon2 = 1e-7;
param_proj.epsilon = epsilon2;
param_proj.A = operatorA;
param_proj.At = operatorAt;
param_proj.y = y;
param_proj.tight = 0;
param_proj.nu = norm(A)^2;
param_proj.verbose = verbose - 1;
f2.prox = @(x,T) proj_b2(x, T, param_proj);
f2.eval = @(x) eps;

% setting the function f1
param_l1.verbose=verbose - 1;
param_l1.tight=1;
f1.prox = @(x, T) prox_l1(x, T, param_l1);
f1.eval = @(x) norm(x,1);   

%% solving the problem

% setting different parameters for the simulation
param_solver.verbose = verbose; % display parameter
param_solver.maxit = 300;       % maximum number of iterations
param_solver.tol = 1e-4;        % tolerance to stop iterating
param_solver.gamma = 1e-2;      % stepsize

% solving the problem
sol = solvep(zeros(N,1) ,{f1, f2}, param_solver);

%% displaying the result
% figure;
plot(1:N, x, 'o', 1:N, sol, 'xr');
legend('Original signal', 'Reconstructed signal');
    

%% Closing the toolbox
close_unlocbox();
