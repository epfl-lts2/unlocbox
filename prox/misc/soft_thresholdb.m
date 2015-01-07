function [sz] = soft_thresholdb(z,T)
%SOFT_THRESHOLDB soft thresholding for mixed sparsity
%   Usage:  soft_threshold(z,T)
%
%   Input parameters:
%         z     : Input signal.
%         T     : Threshold.
%   Output parameters:
%         sz    : Soft thresholded signal.
%

% Nathanael Perraudin
% Date: 14 May 2013

    % handle the size
    size_z=size(z);
    z=z(:);
    T=T(:);
    
    % Precaution on T
    T(T==inf)=realmax;
    

    % This soft thresholding function support complex number
    sz = max(abs(z)-T.*abs(z),0)./(max(abs(z)-T.*abs(z),0)+T.*abs(z)+double(abs(z)<eps)).*z;
    

    % Handle the size
    sz=reshape(sz,size_z);
end
