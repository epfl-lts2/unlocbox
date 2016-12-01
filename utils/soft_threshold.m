function [sz] = soft_threshold(z,T)
%SOFT_THRESHOLD soft thresholding
%   Usage:  sz = soft_threshold(z,T);
%
%   Input parameters:
%         z     : Input signal
%         T     : Threshold
%                 if T is a vector, then thresholding is applied component-wise
%
%   Output parameters:
%         sz    : Soft thresholded signal
%   
%   This function soft thresholds z by T. It can handle complex input z.

% Nathanael Perraudin
% Date: 14 May 2013


size_z = size(z);

if all(T==0)
    sz = z;  %identity
elseif any(T<0)
    error('Threshold value(s) cannot be negative!')
elseif isscalar(T)  %for scalar threshold it is faster to compute it like this
    % handle the size
    z=z(:);
    
    % This soft thresholding function only supports real signal
    % sz = sign(z).*max(abs(z)-T, 0);
    
    % This soft thresholding function supports complex numbers
    sz = max(abs(z)-T,0)./(max(abs(z)-T,0)+T).*z;
    
else  %for vector threshold(s) it is faster to compute it like this
    % handle the size
    z=z(:);
    T=T(:);
    
    % This soft thresholding function supports complex numbers
    % sz = max(abs(z)-T,0)./(max(abs(z)-T,0)+T).*z;
    aux = max(abs(z)-T,0);
    sz = aux./(aux+T).*z;

end

% Handle the size
sz = reshape(sz,size_z);
