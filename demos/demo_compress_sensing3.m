%DEMO_COMPRESS_SENSING3 Compress sensing example using grouped L12 norm
%   
%   We present a compress sensing example solved with the douglas rachford
%   solver. The particularity of this example is the use of a mixed norm. We
%   do not only know the the signal is sparse, we also know that the
%   sparse coefficients are grouped.
%
%   The problem can be expressed as this
%
%   ..   argmin || x||_{2,1} s.t ||b-Ax||_2 < epsilon
%
%   .. math:: arg \min_x \| x\|_{2,1} \hspace{1cm} such \hspace{0.25cm}  that \hspace{1cm} \|b-Ax\|_2 \leq \epsilon
%  
%   Where b are the measurements and A the measurement matrix.
%
%   We set 
%
%   * $f_1(x)=||x||_{2,1}$
%     We define the prox of $f_1$ as: 
%
%     .. prox_{f1,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma ||z||_{2,1}
%
%     .. math:: prox_{f1,\gamma} (z) = arg \min_{x} \frac{1}{2} \|x-z\|_2^2  +  \gamma \|z\|_{2,1}
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
%   The theoretical number of measurements $M$ is computed with respect of the size of the
%   signal $N$ and the sparsity level $K$:
%   
%   .. M=K*max(4,ceil(log(N)))   
%
%   .. math:: M=K \max\left( 4 , \text{ceil}(\log(N))\right).
%
%   Since we add some new information, we will try to reduce the number of
%   measurements by a factor p:
%
%   .. M=K*max(4/p,ceil(log(N)/p))   
%
%   .. math:: M=K \max\left( \frac{4}{p} , \text{ceil}(\frac{\log(N)}{p})\right).
%
%   With this number of measurements, we hope that the algorithm will
%   perform  a perfect reconstruction.
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
%      measurements is M=900, the length of the signal N=5000, K=100, p=4. This is
%      equivalent to a compression ratio of 16.67. The elements are grouped
%      by 10.
%
%   References: combettes2007douglas bach2011optimization combettes2011proximal

 
% Author: Nathanael Perraudin
% Date: november 2012


%% Initialisation

clear;
close all;

% Loading toolbox
init_unlocbox();

verbose = 2; % verbosity level

%% Creation of the problem

g = 10;         % number of element per group
N = 5000;       % Size of the signal
K = g * 10;     % Sparsity level
p = 4;          % Gain with respect to the "traditional compression ratio"
R = max(ceil(4/p),ceil(log(N)/p));    % Constant 
fprintf('The compression ratio is: %g\n',N/(R*K));

% Mesurements matrix
A = randn(R * K, N);

% Create a K sparse signal
x = zeros(N, 1);
I2 = randperm(N/g)*g;
I = zeros(size(x));

for i=0:N/g-1
    I(i*g+1:(i+1)*(g)) = I2(i+1) * ones(g,1) - (1:g)' + 1;
end

x(I(1:K)) = randn(K, 1); % take a normal distribution (adapted to L12 norm)
x = x / norm(x);

% Create groups;
g_d = 1:N;
g_t = g*ones(1, N/g);

% Measurements
y = A * x;

%% Defining proximal operators

% Define the prox of f2 see the function proj_B2 for more help
operatorA = @(x) A * x;
operatorAt = @(x) A' * x;
epsilon2 = 1e-5;
param_proj.epsilon = epsilon2;
param_proj.A = operatorA;
param_proj.At = operatorAt;
param_proj.y = y;
param_proj.tight = 0;
param_proj.nu = norm(A)^2;
param_proj.verbose = verbose - 1;
f2.prox = @(x,T) proj_b2(x, T, param_proj);
f2.eval = @(x) norm(A*x - y)^2;

% setting the function f1
param_l21.verbose = verbose - 1;
param_l21.g_d = g_d;
param_l21.g_t = g_t;
f1.prox = @(x, T) prox_l21(x, T, param_l21);
f1.eval = @(x) norm_l21(x,g_d,g_t);   

%% solving the problem

% setting different parameters for the simulation
param_solver.verbose = verbose; % display parameter
param_solver.maxit = 300;       % maximum number of iterations
param_solver.tol = 1e-4;        % tolerance to stop iterating
param_solver.gamma = 1e-2;      % stepsize

% solving the problem
sol = solvep(zeros(N,1), {f1, f2}, param_solver);

%% displaying the result
figure;
plot(1:N, x, 'o', 1:N, sol, 'xr');
legend('Original signal', 'Reconstructed signal');
    
%% Closing the toolbox
close_unlocbox();
