function y = tv_norm3d(u,wx, wy, wz)
%TV_NORM3D 3 Dimentional TV norm
%   Usage:  y = tv_norm3d(x)
%           y = tv_norm3d(x, wx, wy, wz )
%
%   Input parameters:
%         x     : Input data (3 dimentional matrix)
%         wx    : Weights along x
%         wy    : Weights along y
%         wz    : Weights along z
%
%   Output parameters:
%         y   : Norm
%
%   Compute the 3-dimentional TV norm of x. If the input I is a 4
%   dimentional signal. This function will compute the norm of all cubes
%   and return a vector of norms.  

% Author: Nathanael Perraudin
% Date:   1 February 2014


warning('This function will be deleted in future release of matlab, use norm_tv3d instead')

if nargin>1
    [dx, dy, dz] = gradient_op3d(u,wx, wy, wz);
else
    [dx, dy, dz] = gradient_op3d(u);
end
    
temp = sqrt(abs(dx).^2 + abs(dy).^2 + abs(dz).^2);
% y = sum(temp(:));

% This allows to return a vector of norms
y = reshape(sum(sum(sum(temp,1),2),3),[],1);


end
