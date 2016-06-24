function [dx, dy] = gradient_op(I, wx, wy)
%GRADIENT_OP 2 Dimensional gradient operator
%   Usage:  [dx, dy] = gradient_op(I)
%           [dx, dy] = gradient_op(I, wx, wy)
%
%   Input parameters:
%         I     : Input data 
%         wx    : Weights along x
%         wy    : Weights along y
%
%   Output parameters:
%         dx    : Gradient along x
%         dy    : Gradient along y
%
%   Compute the 2-dimensional gradient of I. If the input I is a cube. This
%   function will compute the gradient of all image and return two cubes.
%
%   See also: gradient_op3d gradient_op1d div_op laplacian_op

% Author: Nathanael Perraudin
% Date:   1 February 2014

dx = [I(2:end, :,:)-I(1:end-1, :,:) ; zeros(1, size(I, 2),size(I, 3))];
dy = [I(:, 2:end,:)-I(:, 1:end-1,:) , zeros(size(I, 1), 1,size(I, 3))];

if nargin>1
    dx = dx .* wx;
    dy = dy .* wy;
end

end
