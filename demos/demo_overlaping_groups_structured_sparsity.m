%DEMO_OVERLAPING_GROUPS_STRUCTURED_SPARSITY Demonstration for the use of overlaping structured sparsity
%
%   Here we try to deblur an image through a deconvolution problem. The
%   convolution operator is the blur. 
%   In order to get better results, we try to group some pixel of the
%   wavelet transform.
%
%   The group of 4 pixels are made in H(x) like this:
%
%   For a 4 by 4 image:: 
%
%       1122        1221
%       1122        3443
%       3344        3443
%       3344        1221
%
%   This is probably not the best way to do it.
%
%   The problem can be expressed as this
%
%   ..   argmin  ||Ax-b||^2 + tau*||H(x)||_12
%
%   .. math:: arg \min_x \|Ax-b\|^2 + \tau \|H(x)\|_{12}
%  
%   Where b is the degraded image, I the identity and A an operator representing the mask.
%
%   H is a linear operator projecting the signal in a sparse
%   representation. Here we worked with wavelet. 
%
%   Warning! Note that this demo require the rwt(RICE WAVELET TOOLBOX) to work.
%
%   We set 
%
%   * $f_1(x)=||H(x)||_{12}$
%     We define the prox of $f_1$ as: 
%
%     .. prox_{f1,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma  ||H(z)||_12
%
%     .. math:: prox_{f1,\gamma} (z) = arg \min_{x} \frac{1}{2} \|x-z\|_2^2  +  \gamma \|H(z)\|_{12}
%
%   * $f_2(x)=||Ax-b||_2^2$
%     We define the gradient as: 
%
%     .. grad_f(x) = 2 * A^*(Ax-b)
%
%     .. math:: \nabla_f(x) = 2 A^*(x-b)
%
%   Results
%   -------
%
%   .. figure::
%
%      Original image
%
%      This figure shows the original cameraman image. 
%
%   .. figure::
%
%      Depleted image
%
%      This figure shows the image after the application of the blur.
%
%   .. figure::
%
%      Reconstruted image
%
%      This figure shows the reconstructed image thanks to the algorithm.
%
%   References: raguet2011generalized


% Author: Nathanael Perraudin
% Date: sept 30 2011
%


%% Initialisation

clear all;
close all;

% Loading toolbox
init_unlocbox();

verbose = 2;    % verbosity level

%% Defining the problem

tau = 0.001; %parameter for the problem

% Original image
im_original=barbara();


% Create a blur
sigma = 0.05;
[x, y] = meshgrid(linspace(-1, 1, length(im_original)));
r = x.^2 + y.^2;
G = exp(-r/(2*sigma^2));

A=@(x) real(ifft2(fftshift(G).*(fft2(x))));
At=@(x) real(ifft2(conj(fftshift(G)).*(fft2(x))));

% Blur the original image
b = A(im_original);

% The group of 4 pixels are made in H(x) like this
%
%   For a 4 by 4 image 
%   1122        1221
%   1122        3443
%   3344        3443
%   3344        1221

% Create the group (this is a bad code)
    % -------------------------------------------- %
    g_d1=zeros(size(b,1)*size(b,2),1);
    indice=1:length(g_d1);
    mat_indice=reshape(indice,size(b))';


    k=size(b,1)/2;
    l=size(b,2)/2;
    for i=1:k
       for j=1:l
           g_d1((i-1)*4*l+4*(j-1)+1:(i-1)*4*l+4*(j-1)+4) = ...
            reshape(mat_indice(2*(i-1)+1:2*(i-1)+2,2*(j-1)+1:2*(j-1)+2),4,1);                                                                                                                                           
       end
    end

    temp=reshape(g_d1,size(b));
    temp=[temp(2:end,:);temp(1,:)];
    temp=[temp(:,2:end),temp(:,1)];
    g_d2=zeros(size(b,1)*size(b,2),1);
    for i=1:k
       for j=1:l
           g_d2((i-1)*4*l+4*(j-1)+1:(i-1)*4*l+4*(j-1)+4) = ...
            reshape(temp(2*(i-1)+1:2*(i-1)+2,2*(j-1)+1:2*(j-1)+2),4,1);                                                                                                                                           
       end
    end


    g_t1=4*ones(length(g_d1)/4,1);
    g_t2=g_t1;
    
    % -------------------------------------------- %
  
param_l12.g_t=[g_t1'; g_t2'];
param_l12.g_d=[g_d1'; g_d2'];

% Wavelet operator for the l12 norm
L=6;
% h = daubcqf(2);
% A2 = @(x) mdwt(x,h,L);
% A2t = @(x) midwt(x,h,L);
A2 = @(x) fwt2(x,'db2',L);
A2t = @(x) ifwt2(x,'db2',L);

param_l12.verbose=1;

f.prox=@(x, T) x+ A2t(prox_l12(A2(x), T*tau, param_l12)-A2(x));
f.eval=@(x) tau*norm_l12(A2(x));   


%% solving the problem

% setting different parameter for the simulation
param_solver.verbose = verbose;     % display parameter
param_solver.maxit = 50;            % maximum iteration
param_solver.tol = 10e-7;           % tolerance to stop iterating
param_solver.gamma = 0.5;           % stepsize (beta is equal to 2)
param_solver.method = 'FISTA';      % desired method for solving the problem

% solving the problem
sol=rlr(b,f,A,At,param_solver);


%% displaying the result
imagesc_gray(im_original, 1, 'Original image');
imagesc_gray(b, 2, 'Depleted image');
imagesc_gray(sol, 3, 'Reconstructed image');
    

%% Closing the toolbox
close_unlocbox();