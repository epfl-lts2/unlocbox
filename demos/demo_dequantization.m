%DEMO_DEQUANTIZATION Dequantization demo
%   This demo shows how a quantized signal, sparse in the DCT domain, can be dequantized
%   solving a convex problem using Douglas-Rachford algorithm
%
%   Suppose signal y has been quantized. In this demo we use quantization levels 
%   that are uniformly spread between the min. and max. value of the
%   signal. The resulting signal is y_Q.
%
%   The problem can be expressed as 
%
%   ..   argmin_x  || x ||_1     s.t.   ||Dx - y_Q||_infty <= alpha/2 
%
%   .. math:: arg \min_x  \| x \|_{1} \text{ s.t. } \|Dx - y_Q\|_\infty \leq \frac{\alpha}{2} 
%
%   where D is the synthesis dictionary (DCT in our case) and $\alpha$ is the
%   distance between quantization levels. The constraint basically
%   represents the fact that the reconstructed signal samples must stay
%   within the corresponding quantization stripes.
%
%   After sparse coordinates are found, the dequantized signal
%   is obtained simply by synthesis with the dictionary.   
%
%   The program is solved using Douglas-Rachford algorithm. We set 
%
%   * $f_1(x)=||x||_{1}$. Its respective prox is the soft thresholding operator.
%
%   * $f_2(x)=i_C$ is the indicator function of the set C, defined as
%
%     .. C = { x | ||Dx - y_Q||_infty <= alpha/2 } 
%
%     .. math:: C = \{ x | \|Dx - y_Q\|_\infty <= \frac{\alpha}{2} \}
%
%   Its prox is the orthogonal projection onto that set, which is realized
%   by entry-wise 1D projections onto the quantization stripes. This is
%   realized for all the entries at once by function proj_box.
%
%   As an alternative, setting algorithm = 'LP' switches to computing the
%   result via linear programming (requires Matlab optimization toolbox).
%
%   Results
%   -------
%
%   .. figure::
%
%      Original, quantized and dequantized signals
%
%       
%
%   .. figure::
%
%      Quantization error and error of reconstruction (i.e. original - reconstr.)
%
%      
%
%   .. figure::
%
%      Coefficients of original and reconstructed signals
%
%      
%
%   References: combettes2007douglas

% Authors: Pavel Rajmic, Pavel Záviška
% Date: August 2015



clc
clear
close all

%% Parameters

N = 64; %length of the signal
k = 5; %sparsity of the signal in the dictionary
d = 10; %number of quantization levels

% 'DR' ... Douglas-Rachford algorithm
% 'LP' ... linear programming (requires Matlab optimization toolbox)
algorithm =  'DR' %'LP' %
                    


%% Dictionary and the original sparse signal
% A = dctmtx(N)';  % dictionary in which the signal will be sparse = inverse DCT
A = @(x) idct(x);  % dictionary in which the signal will be sparse = inverse DCT
At = @(x) dct(x); %its inverse (= transpose)

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
% y = A(:, support)*x_support;
y = A(x);

%or just load a fixed dataset
% clear
% load('dequantization_dataset_01')

% show coefficients
fig_coef = figure;
h_coefs_orig = bar(x);
hold on;
title('Original coefficients of sparse signal');



%% Quantize signal
%Determine the max and min value of the signal. The quantization
%levels are spread uniformly between them
min_y = min(min(y));
max_y = max(max(y));
range = max_y - min_y;

quant_step   = range / (d-1); %quantization step
dec_bounds = (min_y+quant_step/2) : quant_step : (max_y-quant_step/2);  %decision boundaries lie in the middle
quant_levels = min_y : quant_step : max_y; %quantization levels

%use Matlab function to quantize signal
[index, y_quant] = quantiz(y, dec_bounds, quant_levels);
y_quant = y_quant';

%constraints for the signal samples;
%quantized signal lies on the quant. levels and the true signal cannot be outside the boundaries
lower_dec_bounds = y_quant - (quant_step/2);
upper_dec_bounds = y_quant + (quant_step/2);

%the highest and lowest boundary coincide precisely with the quantization levels!:
min_quant_level = quant_levels(1);
max_quant_level = quant_levels(end);
upper_dec_bounds(upper_dec_bounds > max_quant_level) = max_quant_level;
lower_dec_bounds(lower_dec_bounds < min_quant_level) = min_quant_level;


%% Show original vs. quantized signal
%define colors
grey = 0.6* ones(1,3);
lightgrey = 0.8* ones(1,3);
black = [0,0,0];
blue = [0.251,0.584,0.808];
orange = [0.851,0.325,0.098];
green = [0 1 0];

%plot of the signals
fig_time = figure;
h_orig = plot(y, '.-', 'Color', blue);
hold on;
h_quant = plot(y_quant, '.-', 'Color', orange);
title('Original and quantized signal');

%plot signal constraints
h_signal_constr = plot(upper_dec_bounds, 'Color', lightgrey);
plot(lower_dec_bounds, 'Color', lightgrey);

for j=1:d %quantization levels
    yPos = quant_levels(j);
    h_q_lev = plot(ones(1,N) * yPos, 'Color', grey);
end

for j=1:(d-1) %decision boundaries
    yPos = dec_bounds(j);
    h_dec_b = plot(ones(1,N) * yPos, ':', 'Color', grey); 
end

axis tight

%legend
uistack(h_orig,'top');
uistack(h_quant,'top');
h_legend_time = legend([h_orig, h_quant, h_q_lev, h_dec_b, h_signal_constr ], 'original', 'quantized', 'quantiz. levels', 'decision bounds', 'signal constraints');

%quanization noise
fig_quant_noise = figure;
quant_noise = y - y_quant;
h_quant_noise = plot(quant_noise, 'Color', blue);
hold on;
title('Quantization noise');


%% Sparse dequantization
switch algorithm
    case 'LP'  %dequantization using linear programming (doubles the number of variables)
        f = (ones([1 2*N]));
        
        b = [-lower_dec_bounds; upper_dec_bounds];
        Amatrix = A(eye(N)); %generate explicit dictionary matrix
        % A_ = [A -A; -A A];
        A_ = [Amatrix -Amatrix; -Amatrix Amatrix];
        lb = zeros(2*N,1);  %all variables must be nonnegative
        
        w = linprog(f,A_,b,[],[],lb);  %l1-minimization via lin. program
        
        uv = reshape(w,N,2);      %split w into two
        u = uv(:, 1);
        v = uv(:, 2);
        x_reconstructed = v - u;    %sparse vector is determined by subtracting the two non-negative
        y_dequant = A(x_reconstructed);    %reconstruct signal
        
        sol = x_reconstructed;
        
    case 'DR'  % Douglas-Rachford
        
        param.lower_lim = lower_dec_bounds;
        param.upper_lim = upper_dec_bounds;
        
        indi_thr = 10e-5; % threshold for identifying the point to lie in the set
        
        f1.eval = @(x) norm(x,1);
        f1.prox = @(x,T) sign(x).*max(abs(x)-T, 0);
        
        f2.eval = @(x) 1 / (1 - ( any((A(x(:))-param.upper_lim)>indi_thr)) || any((A(x(:))-param.lower_lim) < -indi_thr) ) - 1; %zero if x is inside boundaries, Inf otherwise
        f2.prox = @(x,T) At(proj_box(A(x),[],param)); %box projection in the signal space (thanks to DCT being orthogonal)
        
        %%%%%%%% UNLOCBOX version %%%%%%%%%%%%%
        % setting different parameter for the simulation
        paramsolver.verbose = 5;  % display parameter
        paramsolver.maxit = 300;        % maximum iteration
        paramsolver.tol = 10e-7;        % tolerance to stop iterating
        paramsolver.lambda = 1;   % step for DR algorithm (default 1)
        paramsolver.gamma = 1e-2;        % here threshold for soft thresholding
        
        [sol, info] = douglas_rachford(At(y_quant), f1, f2, paramsolver);
        info
        
        sol = f2.prox(sol,[]); %final projection into the constraints
        y_dequant = A(sol);
        
%         pause
%         
%         %%%%%%%% manual version %%%%%%%%%%%%%      
%         %starting point
%         DR_y = A(y_quant);
%         DR_x_old = DR_y;
%         
%         relat_change_coefs = 1;
%         relat_change_obj = 1;
%         cnt = 1; %iteration counter
%         obj_eval = [];
%         
% %         while relat_change_coefs > 0.00001
%         while relat_change_obj > paramsolver.tol
%             % DR: gamma = 1
%             DR_x = f2.prox(DR_y,[]);
%             obj_eval = [obj_eval, f1.eval(DR_x) + f2.eval(DR_x)]; %record values of objective function
%             DR_y = DR_y + paramsolver.lambda*(f1.prox(2*DR_x-DR_y, paramsolver.gamma)-DR_x);
%             if cnt > 1
%                 relat_change_coefs = norm(DR_x-DR_x_old) / norm(DR_x_old);
%                 relat_change_obj = norm(obj_eval(end) - obj_eval(end-1)) / norm(obj_eval(end-1));
%                 if paramsolver.verbose > 1
%                     fprintf('  relative change in coefficients: %e \n', relat_change_coefs);
%                     fprintf('  relative change in objective fun: %e \n', relat_change_obj);
%                     fprintf('\n');
%                 end
%             end
%             DR_x_old = DR_x;
%             cnt = cnt + 1;
%             
%         end
%         
%         DR_x = f2.prox(DR_y); %final projection into the constraints
%         y_dequant = A(DR_x); %dequantized signal
%         
%         disp(['Finished after ' num2str(cnt) ' iterations.'])
%         
%         
%         %compare UNLOCBOX with manual
%         figure
%         plot([y_dequant A(sol)])
%         norm(y_dequant - A(sol))
%         title('UNLOCBOX vs. manual solution')
%         
%         %behaviour of objective through iterations
%         figure
%         plot(obj_eval)
%         title('Objective function value (after projection into constraints in each iteration)')
    
end




%% Show results

figure(fig_time)
h_dequant = plot(y_dequant, '.-', 'Color', green);
uistack(h_dequant,'top');
delete(h_legend_time)
h_legend_time = legend([h_orig, h_quant, h_dequant, h_q_lev, h_dec_b, h_signal_constr ], 'original', 'quantized', 'dequantized', 'quantiz. levels', 'decision bounds', 'signal constraints');
title('Original, quantized and dequantized signals')

%quantization and reconstruction errors
figure(fig_quant_noise)
h_dequant_error = plot(y-y_dequant, 'Color', green);
title('Quantization error and error of reconstruction (i.e. original - reconstr.)');
axis tight
legend([h_quant_noise h_dequant_error], 'Quantizat. error', 'Error of reconstr.')

%coefficients of reconstructed signal
figure(fig_coef)
hold on
h_coefs_dequant = bar(sol,'FaceColor',green);
title('Coefficients of original and reconstructed signals');
legend([h_coefs_orig h_coefs_dequant], 'Coefs of orig. signal', 'Coefs of dequant. signal')
axis tight


figure(fig_time)
axis tight