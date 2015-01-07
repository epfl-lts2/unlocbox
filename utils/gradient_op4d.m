function [dx, dy, dz, dt] = gradient_op4d(I, wx, wy, wz, wt)
%GRADIENT_OP4D 4 Dimentional gradient operator
%   Usage:  [dx, dy, dz, dt] = gradient_op4d(I)
%           [dx, dy, dz, dt] = gradient_op4d(I, wx, wy, wz, wt)
%
%   Input parameters:
%         I     : Input data 
%         wx    : Weights along x
%         wy    : Weights along y
%         wz    : Weights along z
%         wt    : Weights along t
%
%   Output parameters:
%         dx    : Gradient along x
%         dy    : Gradient along y
%         dz    : Gradient along z
%         dt    : Gradient along t
%
%   Compute the 4-dimentional gradient of I. If the input I has 5
%   dimentions. This function will compute the gradient of all 4
%   dimentional cubes and return 4 5-dimentionals signals
%
%   See also: gradient_op gradient_op1d div_op laplacian_op gradient_op3d

% Author: Nathanael Perraudin
% Date:   25 April 2014

dx = [ I(2:end, :, :, :, :) - I(1:end-1, :, :, :, :); ...
       zeros(1, size(I,2), size(I,3), size(I,4), size(I,5)) ];
dy = [ I(:, 2:end, :, :, :) - I(:, 1:end-1, :, :, :) , ...
       zeros(size(I,1), 1, size(I,3),size(I,4),size(I,5))];
dz = cat(3, I(:, :, 2:end, :, :) - I(:, :, 1:end-1, :, :) , ...
         zeros(size(I,1), size(I,2), 1, size(I,4), size(I,5)));
dt = cat(4, I(:, :, :, 2:end, :) - I(:, :, :, 1:end-1, :) , ...
         zeros(size(I,1), size(I,2), size(I,3), 1, size(I,5)));
     
if nargin>1
    dx = dx .* wx;
    dy = dy .* wy;
    dz = dz .* wz;
    dt = dt .* wt;
end

end
