%DEMO_DOUGLAS_RACHFORD  Example of use of the douglas_rachford solver
%
%   We present an example of the douglas_rachford solver through an image
%   reconstruction problem.
%   The problem can be expressed as this
%
%   ..  argmin ||x||_TV s.t ||b-Ax||_2 < epsilon
%
%   .. math:: arg \min_x \|x\|_{TV} \hspace{1cm} such \hspace{0.25cm}  that \hspace{1cm} \|b-Ax\|_2 \leq \epsilon
%  
%   Where b is the degraded image, I the identity and A an operator representing the mask.
%
%   Note that the constraint can be inserted in the objective function
%   thanks to the help of the indicative function. Then we recover the
%   general formulation used for the solver of this toolbox.
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
%      that 85% of the pixels have been removed.
%
%   .. figure::
%
%      Reconstructed image
%
%      This figure shows the reconstructed image thanks to the algorithm.
%   
%   References: combettes2011proximal

% Author: Nathanael Perraudin
% Date: November 2012
%


%% Initialisation

clear;
close all;

init_unlocbox;

verbose = 2;    % Verbosity level

%% Creation of the problem

% Original image
im_original = lena();   

% Creating the problem
A = rand(size(im_original));
A = A > 0.85;

% Depleted image
b = A .* im_original;

%% Defining proximal operators

% Define the prox of f2 see the function proj_B2 for more help
operatorA = @(x) A.*x;
epsilon2 = 0;
param_proj.epsilon = epsilon2;
param_proj.A = operatorA;
param_proj.At = operatorA;
param_proj.y = b;
param_proj.verbose = verbose - 1;
f2.prox=@(x,T) proj_b2(x, T, param_proj);
f2.eval=@(x) eps;

f2.prox = @(x,T) (x - A.*x) + A.*b;


% setting the function f1 (norm TV)
param_tv.verbose = verbose - 1;
param_tv.maxit = 50;

f1.prox = @(x, T) prox_tv(x, T, param_tv);
f1.eval = @(x) norm_tv(x);   

%% solving the problem

% setting different parameter for the simulation
paramsolver.verbose = verbose;  % display parameter
paramsolver.maxit = 100;        % maximum iteration
paramsolver.tol = 10e-7;        % tolerance to stop iterating
paramsolver.gamma = 0.1;        % stepsize

% To see the evolution of the reconstruction
fig = figure(100);
paramsolver.do_sol = @(x) plot_image(x,fig);

sol = douglas_rachford(b,f1,f2,paramsolver);

close(100);

%% displaying the result
imagesc_gray(im_original, 1, 'Original image');
imagesc_gray(b, 2, 'Depleted image');
imagesc_gray(sol, 3, 'Reconstructed image');
    

%% Closing the toolbox
close_unlocbox();

