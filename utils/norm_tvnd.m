function y = norm_tvnd(u,type,weights)
%NORM_TVND N Dimentional TV norm
%   Usage:  norm_tvnd(x,weights)
%
%   Input parameters:
%         x     : Input data (N dimentional matrix)
%         type  : type ('isotropic' or 'anisotropic') (default 'isotropic')
%         weights: Weights
%   Output parameters:
%         sol   : Norm
%
%   Compute the N-dimentional TV norm of x
%
%   See also: norm_tv norm_tv3d


    if nargin < 2
        type = 'isotropic';
    end

    sz = size(u);
    dim = length(sz);
    
    if nargin<3
        weights = ones(dim,1);
    end


    temp = zeros(sz);

    if strcmp(type,'anisotropic')
        for d = 1:dim
            sz_temp = sz;
            sz_temp(d) = 1;
            tv(d).grad = weights(d)*cat(d,diff(u,1,d),zeros(sz_temp));
            temp = abs(tv(d).grad)+ temp;
        end
        y = sum(temp(:));
    elseif strcmp(type,'isotropic')
        for d = 1:dim
            sz_temp = sz;
            sz_temp(d) = 1;
            tv(d).grad = weights(d)*cat(d,diff(u,1,d),zeros(sz_temp));
            temp = abs(tv(d).grad.^2)+ temp;
        end
        temp = sqrt(temp);
        y = sum(temp(:));
    else
        error('NORMTV_ND: unknown type.')
    end
end
