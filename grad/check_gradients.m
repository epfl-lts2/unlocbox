function is_ok = check_gradients(f_eval, grad, size_input, x0, mode, verbose, symmetric)
% is_ok = check_gradients(f_eval, grad, size_input, x0, mode, verbose, symmetric)
%   f:              function to be checked
%   grad:           gradient of f
%   size_input:     if f_eval(X) is valid, then size_input = size(X);
%   x0:             the point at which the gradient is checked. If not
%                   provided, this is a randn matrix or vector. This is
%                   very useful when additional constraints are used for f,
%                   for example when it is only defined for positive
%                   values.
%   mode:           'directive' or 'bruteforce' (the latter is not
%                   implemented yet!
%   verbose:        0: no output, 1: details
%   symmetric:      TO BE DONE
%
%
% code author: Vassilis Kalofolias

%% TODO
if nargin < 4
    x0 = [];
end
if nargin < 5 || isempty(mode)
    mode = 'directive';
end
if nargin < 6
    verbose = 1;
end
if nargin < 7
    symmetric = 0;
end

if symmetric
    if norm(x0-x0','fro')/norm(x0,'fro') > 0.001
        error('give symmetric input x0')
    end
end


if strcmpi(mode, 'bruteforce')
    
    %TODO
    
    
elseif strcmpi(mode, 'directive')
    
    if isempty(x0)
        x0 = randn(size_input);
        if symmetric
            x0 = x0+x0';
        end
    end
           
    
    delta = randn(size_input) / 1e5;
    if symmetric
        %TODO!!
        %delta = (delta + delta')/2;
    end
    
    f0 = f_eval(x0);
    f1 = f_eval(x0 + delta);
    
    df_approx = (f1 - f0) * 1e5;
    if symmetric
        %TODO!!  
        %df = sum(vec(   g(x0)   .*delta)) * 1e5;
    else
        df = sum(vec(grad(x0).*delta)) * 1e5;
    end
    
    is_ok = abs(df - df_approx) / abs(df) < 1e-3;
end

if verbose
    fprintf('Theoretical derivative: %f, approximation: %f\n', df, df_approx);
    if is_ok
        fprintf('derivative seems OK\n'); 
    else
        fprintf('derivative NOT OK !!\n'); 
    end
end

