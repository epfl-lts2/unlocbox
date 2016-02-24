%DEMO_SOUND_RECONSTRUCTION Sound time in painting demonstration
%
%   Here we solve a sound in-painting problem. The problem can be
%   expressed as this 
%
%   ..   argmin_x  ||A G^* x-b||^2 + tau * || x ||_1
%
%   .. math:: arg \min_x \|A G^*  x-b\|^2 + \tau \| x \|_{1}
%  
%   where $b$ is the signal at the non clipped part,  $A$ an operator
%   representing the mask selecting the non clipped part of the signal and
%   $G^*$ is the Gabor synthesis operation
%
%   Here the general assumption is that the signal is sparse in the Gabor
%   domain!
%   The noiseless particular case of this problem can be epressed as 
%
%   ..   argmin_x  || x ||_1     s. t.   A G^* x = b
%
%   .. math:: arg \min_x  \| x \|_{1} \text{ s.t. } A G^*  x=b
%
%   Warning! Note that this demo requires the LTFAT toolbox to work.
%
%   We set 
%
%   * $f_1(x)=||x||_{1}$
%     We define the prox of $f_1$ as: 
%
%     .. prox_{f1,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma  ||z||_1
%
%     .. math:: prox_{f1,\gamma} (z) = arg \min_{x} \frac{1}{2} \|x-z\|_2^2  +  \gamma \|z\|_{1}
%
%   * $f_2(x)=||Ax-b||_2^2$
%     We define the gradient as: 
%
%     .. grad_f(x) = 2 * G A^*( A G^*x - b )
%
%     .. math:: \nabla_f(x) = 2 * G A^*(AG^*x-b)
%
%   Results
%   -------
%
%   .. figure::
%
%      Original spectrogram
%
%      This figure shows the original spectrogram.
%
%   .. figure::
%
%      Spectrogram of the depleted sound
%
%      This figure shows the spectrogram after the loss of the sample (We loos 75% of the samples.)
%
%   .. figure::
%
%      Spectrogram of the reconstructed sound
%
%      This figure shows the spectrogram of the reconstructed sound thanks to the algorithm.
%
%   References: combettes2007douglas


% Author: Nathanael Perraudin
% Date: sept 30 2011
%


%% Initialisation

clear;
close all;

% Loading toolbox
init_unlocbox();
ltfatstart(); % start the ltfat toolbox

verbose = 2;    % verbosity level
writefile=0;    % writting wav sound

%% Defining the problem

% Original sound
[sound_original, Fs]=gspi();
sound_original = sound_original(1:2^18);



%%
length_sig=length(sound_original); % Put a small number here if you want to proceed only a part a of the signal
sound_part=sound_original(1:length_sig);

% In oder to write the depleted sound somewhere
if writefile
    wavwrite(sound_part,Fs,'original.wav');
end

Mask = rand(size(sound_part))>0.3;

% Depleted sound
sound_depleted = Mask.*sound_part;
sound_depleted(logical(1-Mask)) = randn(sum(1-Mask(:)),1)*mean(abs(sound_part(:)))/5;
if writefile
    wavwrite(sound_depleted,Fs,'depleted.wav');
end

%% Setting proximal operators

tau = 1e-2; % regularization parameter for the problem

% select a gabor frame for a real signal with a Gaussian window
a=64; % size of the shift in time
M=256;% number of frequencies
F=frametight(frame('dgtreal','gauss',a,M));
% Get the framebounds
GB = M/a;

% Define the Frame operators
Psi = @(x) frana(F,x);
Psit = @(x) frsyn(F,x);


% setting the function f2 (l2 norm)
% f2.grad = @(x) 2*Psi(Mask.*(Mask.*(Psit(x)-sound_depleted)));
% f2.eval = @(x) norm(Mask.*Psit(x)-sound_depleted,'fro')^2;
% f2.beta = 2*GB^2;

% noiseless case
f2.prox = @(x,T) Psi( Psit(x) .* (1-  Mask )+ Mask.* sound_depleted );
f2.eval = @(x) eps;


% setting the function f1 (l1 norm of the Gabor transform)
param_l1.verbose = verbose - 1;

f1.prox=@(x, T) prox_l1(x, T*tau, param_l1);
f1.eval=@(x) tau*norm(x,1);   



%% solving the problem


% setting different parameters  for the simulation
param.verbose = verbose; % display parameter
param.maxit = 30; % maximum iteration
param.tol = 10e-5; % tolerance to stop iterating
%param.do_ts = @(x) log_decreasing_ts(x, 10, 0.1, 80);

% Change the stopping criterion to avoid computing the objective function
% every iteration.
param.stopping_criterion = 'rel_norm_primal'; 

sol = Psit(solvep(Psi(sound_depleted),{f1,f2},param));




%% Evaluate the result
snr_in = snr(sound_part,sound_depleted);
snr_fin = snr(sound_part,sol);


fprintf('The SNR of the initial signal is %g dB \n',snr_in);
fprintf('The SNR of the recovered (FB) signal is %g dB \n',snr_fin);



% In order to write the restored sound somewhere
if writefile
    wavwrite(sol,Fs,'restored.wav');
end
%%
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
