%DEMO_LRJS  Example of the lrjs solver
%
%   We present an example of low rank joint sparsity solver. 
%   The problem can be expressed as this
%
%   ..  argmin ||x||_* + ||x||_12 s.t ||b-Ax||_2 < epsilon
%
%   .. math:: arg \min_x \|x\|_{*} + \|x\|_{12} \hspace{1cm} such \hspace{0.25cm}  that \hspace{1cm} \|b-Ax\|_2 \leq \epsilon
%  
%   Where b are the measurments and A the measurment matrix.
%
%   Note that the constraint can be inserted in the objective function
%   thanks to the help of the indicative function. Then we recover the
%   general formulation used for the solver of this toolbox.
%
%
%   Results
%   -------
%
%   .. figure::
%
%      Original matrix
%
%      This is a random joint-sparse and low matrix 
%
%   .. figure::
%
%      Revover matrix
%
%      This figure the recovered matrix thanks to measurements and the
%      algorithms
%  
%   References: golbabaee2012compressed

% Author: Mohammad Golbabaee, Nathanael Perraudin,
% Date:   21 jan. 2012
%


%% Initialisation

clear;
close all;

% Loading toolbox
init_unlocbox();

verbose = 2;    % verbosity level

%% Defining the problem


% parameter of the problems
n1 = 20; % number of lines
n2 = 20; % number of columns
r = 2;   % rank of data matrix
k1 = 10; % Joint-sparsity level
m =  4*r*(n1+n2-r) ; % nb_measurements = 4 * (data degrees of freedom) 


%--- Random measurement matrix--------
A_mtx = 1/sqrt(m) * randn(m, n1*n2);
A = @(x) A_mtx * x;
At = @(x) A_mtx' * x;

% joint-sparse lowrank data generation
x = zeros(n1,n2);
nz_ind = randsample (n1, k1);
x(nz_ind,:) = randn(k1,r)*randn(r,n2);

%-- Compressed sampling-----
y = A(x(:));
z = 0;         %1e-8*randn(nb_meas*J,1); % Additive white gaussian noise.
y = y + z;

%% Parameter for the solver

% general parameters
param_solver.alpha =  sqrt(2*r/k1); % Regularization parameter
epsilon = 0;                        % For noisy data (Fidelity bound)

param_solver.tight_b2 = 0;          % (0 if A is not tight-frame) (1 if A is a tight-frame)
param_solver.nu_b2 = norm(A_mtx)^2; % Upper bound on the spectral norm of the forward operator 
param_solver.verbose = 2;           % Graphical feedback at each main iteration.
param_solver.maxit=200;

% groups
param_solver.g_d=(reshape(reshape(1:n1*n2,n1,n2)',n1*n2,1))';
param_solver.g_t=(n1*ones(n2,1))';

%% Solve the problem

% starting point
x0=zeros(n1,n2);

% solve the problem
xhat=solve_lrjs(x0,y, epsilon, A, At, param_solver);

%% Display the solution

% compute the error
fprintf('The SNR of the solution is %g dB \n',snr(x,xhat));

% display the results
figure(1)
imagesc(x);
title('Original LRJS matrix')
figure(2);
imagesc(xhat);
title('Recovered matrix')
