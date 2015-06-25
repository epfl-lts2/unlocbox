% This demo shows how a quantized signal, sparse in the DCT domain, can be dequantized
% solving a convex problem using Douglas-Rachford algorithm

% 2015 Pavel Rajmic, Pavel Záviška

clc
clear
close all

%% Set parameters

N = 64; %length of the signal
k = 5; %sparsity of the signal in the dictionary
d = 10; %number of quantization levels

% 'DR' ... Douglas-Rachford algorithm
% 'LP' ... linear programming (requires Matlab optimization toolbox)
algorithm =  'DR' %'LP' %
                    


%% Dictionary and the original sparse signal
A = dctmtx(N)';  % dictionary in which the signal will be sparse = inverse DCT
% A = @(x) idct(x);  % dictionary in which the signal will be sparse = inverse DCT

if k>N
    error('Sparsity cannot be greater then signal length/number of atoms');
end

% generate random support
support = randperm(N);
support = support(1:k);

% Determine the coefficient values
x_support =  1 + 3*rand([k 1]);
x_support =  x_support .* (((rand([k 1])>.5)*2)-1); %randomizing signs
x = zeros(N,1);
x(support) = x_support; %complete vector including zero coefs

% synthesize signal
y = A(:, support)*x_support;

% show coefficients
fig_coef = figure;
bar(x);
title('Original coefficients of sparse signal');


%% Quantize signal
%Determine the max and min value of the signal. The quantization
%levels are spread uniformly between them
min_y = min(min(y));
max_y = max(max(y));
range = max_y - min_y;

quant_step = range / (d-1); %quantization step
dec_bounds = (min_y+quant_step/2) : quant_step : (max_y-quant_step/2);  %decision boundaries lie in the middle
quant_levels = min_y : quant_step : max_y; %quantization levels

min_quant_level = quant_levels(1);
max_quant_level = quant_levels(end);

%use Matlab function to quantize signal
[index, y_quant] = quantiz(y, dec_bounds, quant_levels);
y_quant = y_quant';


%% Show original vs. quantized signal
%define colors
grey = [0.6,0.6,0.6];
black = [0,0,0];
blue = [0.251,0.584,0.808];
orange = [0.851,0.325,0.098];
green = [0 1 0];

%plot of the signals
figure
h1 = plot(y, '--', 'Color', blue);
hold on;
h2 = plot(y_quant, 'Color', orange);
title('Original and quantized signal');

for j=1:d %quantization levels
    yPos = quant_levels(j);
    h3 = plot(get(gca,'xlim'), [yPos yPos], 'Color', grey);
end

for j=1:(d-1) %decision boundaries
    yPos = dec_bounds(j);
    h4 = plot(get(gca,'xlim'), [yPos yPos], '--', 'Color', grey); 
end

axis tight

%legend
uistack(h1,'top');
uistack(h2,'top');
legend([h1,h2,h3,h4], 'original', 'quantized', 'quantiz. levels', 'decision bounds');

%quanization noise
figure
quant_noise = y - y_quant;
hr = plot(quant_noise, 'Color', blue);
hold on;
title('Quantization noise');


%% Sparse dequantization
%constraints for the signal samples
lower_dec_bounds = y_quant - (quant_step/2);
upper_dec_bounds = y_quant + (quant_step/2);

%the highest and lowest boundary coincide precisely with the quantization levels!:
upper_dec_bounds(upper_dec_bounds > max_quant_level) = max_quant_level;
lower_dec_bounds(lower_dec_bounds < min_quant_level) = min_quant_level;

switch algorithm
    case 'LP'  %dequantization using linear programming (doubles the number of variables)
        f = (ones([1 2*N]));
        
        b = [-lower_dec_bounds; upper_dec_bounds];
        A_ = [A -A; -A A];
        lb = zeros(2*N,1);  %all variables must be nonnegative
        
        w = linprog(f,A_,b,[],[],lb);  %l1-minimization via lin. program
        
        uv = reshape(w,N,2);      %split w into two
        u = uv(:, 1);
        v = uv(:, 2);
        x_reconstructed = v - u;    %sparse vector is determined by subtracting the two non-negative
        z = A * x_reconstructed;    %reconstruct signal
        
    case 'DR'  % Douglas-Rachford
        
        param.lower_lim = lower_dec_bounds;
        param.upper_lim = upper_dec_bounds;
        
        f1.eval = @(x) norm(x,1);
        f1.prox = @(x,T) sign(x).*max(abs(x)-T, 0);
        
        f2.eval = @(x) 1 / (1 - ( any(dct(x(:))>param.upper_lim) | any(dct(x(:))<param.lower_lim) )) - 1; %zero if x is inside boundaries, Inf otherwise
        f2.prox = @(x,T) dct(proj_box(idct(x),[],param)); %box projection in the signal space (thanks to DCT being orthogonal)
        
        % setting different parameter for the simulation
        paramsolver.verbose = 5;  % display parameter
        paramsolver.maxit = 100;        % maximum iteration
        paramsolver.tol = 10e-7;        % tolerance to stop iterating
        % paramsolver.gamma = 0.1;        % stepsize
        
 %         [sol,info,objective] = douglas_rachford(dct(y_quant), f1, f2, paramsolver);
        
        % DR: lambda
        lambda = 1;
        
        %starting point
        DR_y = dct(y_quant);
        DR_x_old = DR_y;
        
        relat_change = 1;
        cnt = 1; %iteration counter
        
        while relat_change > 0.0001
            % DR: gamma = 1
            DR_x = f2.prox(DR_y,[]);
            DR_y = DR_y + lambda*(f1.prox(2*DR_x-DR_y, 1)-DR_x);
            if cnt > 1
                relat_change = norm(DR_x-DR_x_old) / norm(DR_x_old);
            end
            fprintf('  relative change in coefficients: %e \n', relat_change);
            DR_x_old = DR_x;
            cnt = cnt + 1;
            
        end
        
        DR_x = f2.prox(DR_y); %final projection into the constraints
        z = idct(DR_x); %dequantized signal
        
        disp(['Finished after ' num2str(cnt) ' iterations.'])
        
end


%% Show results
fig_compar = figure;
%time plots of signals
h1 = plot(y, '--', 'Color', blue);
hold on;
h2 = plot(y_quant, 'Color', orange);
h3 = plot(z, 'Color', green);
title('Original, quantized and reconstructed signals');
hold on;

for j=1:d %quantization levels
    yPos = quant_levels(j);
    h4 = plot(get(gca,'xlim'), [yPos yPos], 'Color', grey);
end

for j=1:(d-1) %decision boundaries
    yPos = dec_bounds(j);
    h5 = plot(get(gca,'xlim'), [yPos yPos], '--', 'Color', grey); 
end

%legend
uistack(h1,'top');
uistack(h2,'top');
uistack(h3,'top');
legend([h1,h2,h3,h4,h5], 'original', 'quantized','reconstructed', 'quant. levels', 'decis. bounds.');

axis tight

%coefficients of recontructed signal
figure(fig_coef)
hold on
bar(A'*z,'FaceColor',green);
title('Coefficients of original and reconstructed signals');
axis tight

figure(fig_compar)
