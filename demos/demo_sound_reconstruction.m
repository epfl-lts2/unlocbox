%DEMO_SOUND_RECONSTRUCTION Sound time in painting demonstration
%
%   Here we try to recover missing sample of an sound. The problem can be expressed as this
%
%   ..   argmin  ||Ax-b||^2 + tau*||G(x)||_1
%
%   .. math:: arg \min_x \|Ax-b\|^2 + \tau \|G(x)\|_{1}
%  
%   Where b is the degraded image, I the identity and A an operator representing the mask.
%
%   G is a linear operator projecting the signal in a sparse
%   representation. Here we are working with a Gabor transform. 
%
%   Warning! Note that this demo requires the LTFAT toolbox to work.
%
%   We set 
%
%   * $f_1(x)=||G(x)||_{1}$
%     We define the prox of $f_1$ as: 
%
%     .. prox_{f1,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma  ||G(z)||_1
%
%     .. math:: prox_{f1,\gamma} (z) = arg \min_{x} \frac{1}{2} \|x-z\|_2^2  +  \gamma \|G(z)\|_{1}
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
%      This figure shows the original histogram.
%
%   .. figure::
%
%      Depleted image
%
%      This figure shows the histogram after the loss of the sample (We loos 75% of the samples.)
%
%   .. figure::
%
%      Reconstructed image
%
%      This figure shows the histogram of the reconstructed sound thanks to the algorithm.
%   References: combettes2007douglas


% Author: Nathanael Perraudin
% Date: sept 30 2011
%


%% Initialisation

clear all;
close all;

% Loading toolbox
init_unlocbox();
ltfatstart(); % start the ltfat toolbox

verbose = 2;    % verbosity level
writefile=0;    % writting wav sound

%% Defining the problem

% Original sound
[sound_original, Fs]=gspi();
length_sig=length(sound_original); % Put a small number here if you want to proceed only a part a of the signal
sound_part=sound_original(1:length_sig);

% In oder to write the depleted sound somewhere
if writefile
    wavwrite(sound_part,Fs,'original.wav');
end

% Creating the problem
Mask=rand(size(sound_part));
Mask=(Mask>0.66);
% Depleted sound
sound_depleted=Mask.*sound_part;
if writefile
    wavwrite(sound_depleted,Fs,'depleted.wav');
end

%% Setting proximal operators

tau = 1e-2; % regularization parameter for the problem



% setting the function f2 (l2 norm)
f2.grad = @(x) 2*Mask.*(Mask.*x-sound_depleted);
f2.prox = @(x,T) x-Mask.*x+sound_depleted;
f2.eval = @(x) norm(Mask(:).*x(:)-sound_depleted(:))^2;


% setting the function f1 (l1 norm of the Gabor transform)

% select a gabor frame for a real signal with a Gaussian window
a=64; % size of the shift in time
M=256;% number of frequencies
F=frametight(frame('dgtreal','gauss',a,M));

% Get the framebounds
[GA,GB]= framebounds(F);

% Define the Frame operators
Psi = @(x) frana(F,x);
Psit = @(x) frsyn(F,x);

% tight frame constant
param_l1.nu = GB;

% set parameters
param_l1.verbose = verbose - 1;
param_l1.maxit = 10;
param_l1.A = Psi;
param_l1.At = Psit;

% Since the space of Gabor coefficient is bigger than the space of time
% coefficient
param_l1.tight = 0;

f1.prox=@(x, T) prox_l1(x, T*tau, param_l1);
f1.eval=@(x) tau*norm(Psi(x),1);   



%% solving the problem

% setting different parameters  for the simulation
param.verbose = verbose; % display parameter
param.maxit = 50; % maximum iteration
param.tol = 10e-5; % tolerance to stop iterating
param.gamma = 0.5; % stepsize (beta is equal to 2)
param.method = 'FISTA'; % desired method for solving the problem

sol=forward_backward(sound_depleted,f1,f2,param);

% use of another solver
param.gamma=10;
sol2=douglas_rachford(sound_depleted,f1,f2,param);

%% Evaluate the result
snr_in = snr(sound_part,sound_depleted);
snr_fin = snr(sound_part,sol);
snr_fin2 = snr(sound_part,sol2);

fprintf('The SNR of the initial signal is %g dB \n',snr_in);
fprintf('The SNR of the recovered (FB) signal is %g dB \n',snr_fin);
fprintf('The SNR of the recovered (DG) signal is %g dB \n',snr_fin2);


% In order to write the restored sound somewhere
if writefile
    wavwrite(sol,Fs,'restored.wav');
end

dr=90;

figure(1);
plotframe(F,Psi(sound_part),Fs,dr);
title('Gabor transform of the original sound');

figure(2);
plotframe(F,Psi(sound_depleted),Fs,dr);
title('Gabor transform of the depleted sound');

figure(3);
plotframe(F,Psi(sol),Fs,dr);
title('Gabor transform of the reconstructed sound');
