% Toto je dávkový soubor, který zahrnuje vytvoøené funkce pro gererování øídkého signálu,
% kvantizaci vygenerovaného sigálu a následou dekvantizaci.

%Pavel Rajmic, Pavel Záviška

clc
clear
close all

%% Nastavení
    % V tomto bloku probíhá celé nastavení. Lze nastavit délku 1D signálu
    % (N), øídkost (k) a poèet kvantovacích hladin (d). 
    % Další nastavení se týkají použité báze, ze které se bude signál 
    % generovat (DCT nebo DWT). V pøípadì DWT lze dále nastavit typ 
    % waveletu a hloubku dekompozice. 
    % Nakonec lze nastavit požadovaný dekvantizaèní algoritmus

N = 64; %length of the signal
k = 4; %sparsity of the signal in the dictionary
d = 8; %number of quantization levels

% 'DR' ... Douglas-Rachford algorithm
% 'LP' ... linear programming (requires Matlab optimization toolbox)
algorithm = 'DR'; %'LP' 
                    


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

% show  MOVE?? no
fig_coef = figure;
bar(x);
title('Original coefficients of sparse signal');

%% Quantize signal
%We determine the max and min value of the signal and the quantization
%levels are spread uniformly betwen them
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


%% Show & compare original vs. quantized signal
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
%         s = (quant_step/2)-y_quant;   %constraints for the signal samples
%         t = (quant_step/2)+y_quant;
%constraints for the signal samples
lower_dec_bounds = y_quant - (quant_step/2);
upper_dec_bounds = y_quant + (quant_step/2);

%the highest and lowest boundary coincide with the quantization levels!:
upper_dec_bounds(upper_dec_bounds > max_quant_level) = max_quant_level;
lower_dec_bounds(lower_dec_bounds < min_quant_level) = min_quant_level;

switch algorithm
    case 'LP'   %dequantization using linear programming (doubles the number of variables)
        f = (ones([1 2*N]));
        
        b = [-low_constr; upp_constr];
        A_ = [A -A; -A A];
        lb = zeros(2*N,1);  %all variables must be nonnegative
                
        w = linprog(f,A_,b,[],[],lb);   %l1-minimalizace pomocí lineárního programování
                
        uv = reshape(w,N,2);      %split w into two
        u = uv(:, 1);
        v = uv(:, 2);
        x_reconstructed = v - u;    %sparse vector is determined by subtracting the two non-negative
        z = A * x_reconstructed;    %reconstruct signal
        
    case 'DR'  % dekvantizace pomocí Douglas-Rachfordova algoritmu s projekcí v oblasi signálu 
        
        
        param.lower_lim = lower_dec_bounds;
        param.upper_lim = upper_dec_bounds;
            
        f1.eval = @(x) norm(x,1);
        f1.prox = @(x, T) sign(x).*max(abs(x)-T, 0);
        % Definice funkce mìkkého prahování
%         soft = @(z, T) sign(z).*max(abs(z)-T, 0);
        
%         f2.eval = @(x) 1 / (any(dct(x)>param.upper_lim) | any(dct(x)<param.lower_lim)) - 1; %zero if inside boundaries, Inf otherwise
%         param.upper_lim = 65*ones(N,1);
%         param.lower_lim = 0*ones(N,1);
        f2.eval = @(x) 1 / (1 - ( any((x(:))>param.upper_lim) | any((x(:))<param.lower_lim) )) - 1; %zero if inside boundaries, Inf otherwise
        f2.prox = @(x) proj_box(x,[],param);
        
        % setting different parameter for the simulation
        paramsolver.verbose = 5;  % display parameter
        paramsolver.maxit = 100;        % maximum iteration
        paramsolver.tol = 10e-7;        % tolerance to stop iterating
        % paramsolver.gamma = 0.1;        % stepsize
  
        [sol,info,objective] = douglas_rachford(dct(y_quant), f1, f2, paramsolver);
        
        % Stanovení parametru lambda
        lambda = 1;
        
        %startovaci body? neni to divne?
        x_old = dct(y_quant-quant_step);
        y_new = dct(y_quant-quant_step);

        rel_odchylka = 1;

        while rel_odchylka > 0.0005
%             x = proj_signal_dct(y_quant, y_new, quant_step, min_level, max_level);
            % x = proj_box()
            y_new_idct = idct(y_new);
            [sol, info] = proj_box(y_new_idct,[],param);
            
            % Pøevedení zpìt do oblasti souøadnic
            x = dct(sol);


            koef = (2*x-y_new);
            y_new = y_new + lambda*(soft(koef,1)-x);
        
            rel_odchylka = (norm(x-x_old)/norm(x_old));
            x_old = x;
        end

%         x = proj_signal_dct(y_quant, y_new, quant_step, min_level, max_level);
        y_new_idct = idct(y_new);
            [sol, info] = proj_box(y_new_idct,[],param);
            x = dct(sol);
        
        z = idct(x);        
        
end

% z_fixed = fix_dekvantizace(z, y_quant, dec_bounds);    % posunutí signálu od rozhodovacích úrovní


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

% demo_dekvantizace(y, z_fixed, d);             % zobrazení výsledku dekvantizace  

%coefficients of recontructed signal
figure(fig_coef)
hold on
bar(A'*z,'FaceColor',green);
title('Coefficients of original and reconstructed signals');
axis tight

figure(fig_compar)


%%%%%%%%%%%%
