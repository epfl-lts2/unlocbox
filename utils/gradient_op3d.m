function [dx, dy, dz] = gradient_op3d(I, wx, wy, wz)
%GRADIENT_OP3D 3 Dimentional gradient operator
%   Usage:  [dx, dy, dz] = gradient_op3d(I)
%           [dx, dy, dz] = gradient_op3d(I, wx, wy, wz)
%
%   Input parameters:
%         I     : Input data 
%         wx    : Weights along x
%         wy    : Weights along y
%         wz    : Weights along z
%
%   Output parameters:
%         dx    : Gradient along x
%         dy    : Gradient along y
%         dz    : Gradient along z
%
%   Compute the 3-dimentional gradient of I. If the input I has 4
%   dimentions. This function will compute the gradient of all cubes and
%   return 3 4-dimentionals signals
%
%   See also: gradient_op gradient_op1d div_op laplacian_op

% Author: Nathanael Perraudin
% Date:   1 February 2014

dx = [I(2:end, :, :,:)-I(1:end-1, :, :,:) ;...
      zeros(1, size(I, 2), size(I, 3),size(I, 4))];
dy = [I(:, 2:end, :,:)-I(:, 1:end-1, :,:) , ...
    zeros(size(I, 1), 1, size(I, 3),size(I, 4))];
dz = cat(3, I(:, :, 2:end,:)-I(:, :, 1:end-1,:) , ...
    zeros(size(I, 1),size(I, 2), 1,size(I, 4)));

if nargin>1
    dx = dx .* wx;
    dy = dy .* wy;
    dz = dz .* wz;
end

end
