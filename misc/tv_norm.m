function y = tv_norm(I,wx,wy)
%TV_NORM 2 Dimentional TV norm
%   Usage:  y = tv_norm(x)
%
%   Input parameters:
%         I     : Input data 
%         wx    : Weights along x
%         wy    : Weights along y
%   Output parameters:
%         y     : Norm
%
%   Compute the 2-dimentional TV norm of I. If the input I is a cube. This
%   function will compute the norm of all image and return a vector of
%   norms.

% Author: Nathanael Perraudin
% Date:   1 February 2014

warning('This function will be deleted in future release of matlab, use norm_tv instead')

if nargin>1
    [dx, dy] = gradient_op(I,wx, wy);
else
    [dx, dy] = gradient_op(I);
end    
temp = sqrt(abs(dx).^2 + abs(dy).^2);

%y = sum(temp(:));
y = reshape(sum(sum(temp,1),2),[],1);

end
