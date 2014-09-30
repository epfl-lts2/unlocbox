function [x0 f1,f2,param ] = create_cs_problem( )
%CREATE_CS_PROBLEM create a compress sensing problem

% setting different parameter for the simulation
param.verbose=1; % display parameter
param.maxit=300; % maximum iteration
param.tol=1e-4; % tolerance to stop iterating



N=500; % Size of the signal
K = 100; % Sparsity level
R=max(4,ceil(log(N)));    % Constant 



% Mesurements matrix
A = randn(R*K, N);
param.gamma = 1e-2; % stepsize (beta is equal to 2)

% Create a K sparse signal
x = zeros(N, 1); I = randperm(N);
x(I(1:K)) = randn(K, 1);
x = x/norm(x);
% Measurements
y = A*x;

% Define the prox of f2 see the function proj_B2 for more help
operatorA=@(x) A*x;
operatorAt=@(x) A'*x;
epsilon2=1e-3;
param3.epsilon=epsilon2;
param3.A=operatorA;
param3.At=operatorAt;
param3.y=y;
param3.tight=0;
param3.nu=norm(A)^2;
param3.verbose=0;


% setting the function f2 
f2.prox=@(x,T) proj_b2(x,T,param3);
f2.eval=@(x) norm(A*reshape(x,N,1)-y)^2;



% setting the function f1


param2.verbose=0;
param2.maxit=50;

f1.prox=@(x, T) prox_l1(x, T, param2);
f1.eval=@(x) norm(x,1);   

x0=zeros(size(x));

end

