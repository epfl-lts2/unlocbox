
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



%%
length_sig=length(sound_original); % Put a small number here if you want to proceed only a part a of the signal
sound_part=sound_original(1:length_sig);

% In oder to write the depleted sound somewhere
if writefile
    wavwrite(sound_part,Fs,'original.wav');
end

Mask=(sound_part<0.1)+(sound_part>-0.1);

% Depleted sound
sound_depleted=Mask.*sound_part;
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
f2.grad = @(x) 2*Psi(Mask.*(Mask.*Psit(x)-sound_depleted));
f2.eval = @(x) norm(Mask.*Psit(x)-sound_depleted,'fro')^2;


% setting the function f1 (l1 norm of the Gabor transform)





% set parameters
param_l1.verbose = verbose - 1;


% Since the space of Gabor coefficient is bigger than the space of time
% coefficient


f1.prox=@(x, T) prox_l1(x, T*tau, param_l1);
f1.eval=@(x) tau*norm(x,1);   



%% solving the problem


% setting different parameters  for the simulation
param.verbose = verbose; % display parameter
param.maxit = 50; % maximum iteration
param.tol = 10e-5; % tolerance to stop iterating
param.gamma = 0.5/(GB^2); % stepsize (beta is equal to 2)
param.method = 'FISTA'; % desired method for solving the problem

sol=Psit(forward_backward(Psi(sound_depleted),f1,f2,param));



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
