%DEMO_ADMM Example of use of the ADMM solver
%
%   The demo file present an example of the ADMM (alternating direction
%   method of multipliers) solver. Unfortunately, this method is not fully
%   automatic and the user needs to define the functions in a particular
%   way.
%
%   Please read the paper of Boyd "Distributed Optimization and Statistical
%   Learning via the Alternating Direction Method of Multipliers" to be
%   able to understand this demonstration file. 
%
%   ADMM is used to solve problem of the form
%
%   .. sol = argmin f1(x) + f2(y) such that y=Lx
%
%   .. math::  sol = \min_x f_1(y) + f_2(x) \hspace{1cm} s.t. \hspace{1cm}  y=Lx
%  
%   In this demonstration file, we tackle the following problem
%
%   ..  argmin  tau || z - M x ||_2^2 + || L x ||_1
%
%   .. math:: arg \min_x  \tau \|Mx-z\|_2^2 + \|Lx\|_1 
%
%   where $z$ are the measurements, $W$ the discrete wavelet transform, $M$
%   a masking operator and $\tau$ a regularization parameter. Clearly,
%   setting $Lx=y$ allows to recover the general form for ADMM problem.
%   Contrarily to the other solvers of the UNLocBoX the solver require
%   special proximal operators.
%
%   Here $f_1(x) = \tau \|Mx-z\|_2^2$ would normally take the following
%   proximal operator::
%
%   		f1.prox = @(x, t) ( 1 + tau * t * mask ).^(-1) .* ( x + tau * t * mask.*z);
%   		f1.eval = @(x) tau * norm(mask .* x - z)^2; 
%
%   which correspond to the solution of the following problem
%
%     .. prox_{f1,t} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  t || M x - y ||_2^2
%
%     .. math:: prox_{f1,t} (z) = arg \min_{x} \frac{1}{2} \| x - z \|_2^2  +  t \| M x - y \|_2^2
%
%   However, the ADMM algorithm requires to solve a special proximal
%   operator instead:
%
%     .. prox_{f1,t}^L (z) = argmin_{x} 1/2 || L x - z ||_2^2  +  t || M x - y ||_2^2
%
%     .. math:: prox_{f1,t}^L (z) = arg \min_{x} \frac{1}{2} \| L x - z \|_2^2  +  t \| M x - y \|_2^2
%
%   which is define in MATLAB as::
%
%   		f1.proxL = @(x, t) ( 1 + tau * t * mask ).^(-1) .* ( Lt(x) + tau * t * mask.*z);
%   		f1.prox = @(x, t) ( 1 + tau * t * mask ).^(-1) .* ( x + tau * t * mask.*z);
%   		f1.eval = @(x) tau * norm(mask .* x - z)^2; 
%
%   where *Lt* it the adjoint of the $L$ ( here the inverse wavelet
%   transform) Because the wavelet transform is an orthonormal basis. 
%
%   The function $f_2(y) = \| y \|_1$ is defined in MATLAB as::
%
%   		param_l1.verbose = verbose - 1;
%   		f2.prox = @(x, T) prox_l1(x, T, param_l1);
%   		f2.eval = @(x) norm_l1(L(x));
%   		f2.L = L;
%   		f2.Lt = Lt;
%
%   Note the field *f2.L* and f2.Lt that indicate that the real function
%   function is actually $f_2(Ly) =\| Lx \|_1$.
%   
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
%      This figure shows the image after the application of the mask and addition of the noise. Note that 50% of the pixels have been removed.
%
%   .. figure::
%
%      Reconstructed image
%
%      This figure shows the reconstructed image thanks to the algorithm.
%   
%
%   References: combettes2011proximal boyd2011distributed

% Author: Nathanael Perraudin
% Date: Mai 2015 
%


%% Initialisation

clear;
close all;

% Loading toolbox
init_unlocbox();
ltfatstart();

verbose = 2;    % verbosity level

% Regularization parameter: weight of the fielity term
tau = 50;
% Noise level
sigma = 0.1;
% Percent of missing pixels
p = 50;

%% Defining the problem


% Original image
im_original = barbara();

% Depleted image
mask = rand(size(im_original))>p/100;
z = mask .* im_original + sigma * rand(size(im_original));



%% Defining proximal operators

% Define the wavelet operator
L = @(x)  fwt2(x,'db8',6);
Lt = @(x)  ifwt2(x,'db8',6);

% setting the function tau * || Mx - y ||_2^2  
f1.proxL = @(x, T) (1+tau*T*mask).^(-1) .* (Lt(x)+tau*T*mask.*z);
f1.eval = @(x) tau * norm(mask .* x - z)^2;

% setting the function || L x ||_1 using ADMM to move the operator ot of
% the proximal
param_l1.verbose = verbose - 1;
f2.prox = @(x, T) prox_l1(x, T, param_l1);
f2.eval = @(x) norm(L(x),1);
f2.L = L;
f2.Lt = Lt;
f2.norm_L = 1;


%% solving the problem

% setting different parameter for the solver
paramsolver.verbose = verbose;     % display parameter
paramsolver.maxit = 100;           % maximum number of iterations
paramsolver.tol = 1e-3;            % tolerance to stop iterating
paramsolver.gamma = 1;             % stepsize
% Activate debug mode in order to compute the objective function at each
% iteration.
paramsolver.debug_mode = 1; 
fig=figure(100);
paramsolver.do_sol=@(x) plot_image(x,fig);  

sol = admm(z, f1, f2, paramsolver);

%% displaying the result
imagesc_gray(im_original, 1, 'Original image');
imagesc_gray(z, 2, 'Depleted image');
imagesc_gray(sol, 3, 'Reconstructed image');
    

%% Closing the toolbox
close_unlocbox();
