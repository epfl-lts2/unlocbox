function I = div_op3d(dx, dy, dz, wx, wy, wz)
%DIV_OP3D Divergence operator in 3 dimentions
%   Usage:  I = div_op3d(dx, dy, dz)
%           I = div_op3d(dx, dy, dz, wx, wy, wz)
%
%   Input parameters:
%         dx    : Gradient along x
%         dy    : Gradient along y
%         dz    : Gradient along z
%         wx    : Weights along x
%         wy    : Weights along y
%         wz    : Weights along z
%
%   Output parameters:
%         I     : Output image
%
%   Compute the 3-dimentional divergence of a 3D-image. If a 4 dimentional
%   signal is given, it will compute the divergence of all cubes in the
%   4 diementionals signal.  
%
%   Warning this function compute the divergence operator defined as minus
%   the adjoint of the gradient
%
%   ..      div  = - grad'
%
%   .. math:: \text{div} = - \nabla^*
%
%   See also: gradient_op div_op div_op1d laplacian_op

% Author: Nathanael Perraudin
% Date:   1 February 2014

if nargin > 3
    dx = dx .* conj(wx);
    dy = dy .* conj(wy);
    dz = dz .* conj(wz);
end

I = [dx(1, :, :,:) ; dx(2:end-1, :, :,:) - ...
    dx(1:end-2, :, :,:) ; -dx(end-1, :, :,:)];
I = I + [dy(:, 1, :,:) , dy(:, 2:end-1, :,:) - ...
    dy(:, 1:end-2, :,:) , -dy(:, end-1, :,:)];
I = I + cat(3, dz(:, :, 1,:) , dz(:, :, 2:end-1,:) - ...
    dz(:, :, 1:end-2,:) , -dz(:, :, end-1,:));
end
