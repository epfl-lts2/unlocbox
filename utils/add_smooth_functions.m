function fg = add_smooth_functions(Fg)
%ADD_SMOOTH_FUNCTIONS sum of smooth function
%   Usage: fg = add_smooth_functions(Fg);
%    
%   Input parameters:
%       Fg  : cell array of function (cell array of struct)
%   Output parameters:
%       fg  : function (struct)
%
%   This function takes a cell array of smooth functions and transform into
%   a single smooth function. The array can be empty.
%


if numel(Fg)==0
    fg.grad = @(x) 0;
    fg.eval = @(x) 0;
    fg.beta = 0;
elseif numel(Fg) == 1
    fg = Fg{1};
else
    beta = 0;
    for ii = 1:numel(Fg)
        beta = beta + Fg{ii}.beta;
    end

    fg.grad = @(x) sum_grad(Fg,x);
    fg.eval = @(x) sum_eval(Fg,x);
    fg.beta = beta;
end

if ~(numel(fg.beta)==1)
    error('f.beta must be a scalar!')
end

end



function curr_grad = sum_grad(Fg,x)

    curr_grad = 0;
    for ii = 1:numel(Fg)
        curr_grad = curr_grad + Fg{ii}.grad(x);
    end
    
end

function curr_eval = sum_eval(Fg,x)

    curr_eval = 0;
    for ii = 1:numel(Fg)
        curr_eval = curr_eval + Fg{ii}.eval(x);
    end

end