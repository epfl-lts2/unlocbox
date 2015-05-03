function [sz] = soft_threshold(z,T)
%SOFT_THRESHOLD soft thresholding
%   Usage:  sz = soft_threshold(z,T);
%
%   Input parameters:
%         z     : Input signal
%         T     : Threshold
%   Output parameters:
%         sz    : Soft thresholded signal
%   
%   This function soft threshold z by T. It can handle complex numbers.

% Nathanael Perraudin
% Date: 14 May 2013

    % handle the size
    if T>0
        size_z=size(z);
        z=z(:);
        T=T(:);

        % This soft thresholding function only supports real signal
        % sz= sign(z).*max(abs(z)-T, 0);

        % This soft thresholding function supports complex numbers
        sz = max(abs(z)-T,0)./(max(abs(z)-T,0)+T).*z;

        % Handle the size
        sz=reshape(sz,size_z);
    else
        sz = z;
    end
end
