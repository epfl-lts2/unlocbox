function y = tv_normnd(u,weights)
%TV_NORMND N Dimentional TV norm
%   Usage:  tv_normnd(x,weights)
%
%   Input parameters:
%         x     : Input data (N dimentional matrix)
%         weights: Weights
%   Output parameters:
%         sol   : Norm
%
%   Compute the N-dimentional TV norm of x


warning('This function will be deleted in future release of matlab, use norm_ndtv instead')

sz = size(u);
        anisotropic = 1;
        dim = length(sz);
        temp = zeros(sz);
        if nargin<2
            weights = ones(dim,1);
        end
    if anisotropic == 1
        for d = 1:dim
            sz_temp = sz;
            sz_temp(d) = 1;
            tv(d).grad = weights(d)*cat(d,diff(u,1,d),zeros(sz_temp));
            temp = abs(tv(d).grad)+ temp;
        end
        y = sum(temp(:));
    else
        for d = 1:dim
            sz_temp = sz;
            sz_temp(d) = 1;
            tv(d).grad = weights(d)*cat(d,diff(u,1,d),zeros(sz_temp));
            temp = abs(tv(d).grad.^2)+ temp;
        end
        temp = sqrt(temp);
        y = sum(temp(:));
    end
end
