function I = div_op(dx, dy, wx, wy)
%DIV_OP Divergence operator in 2 dimensions
%   Usage:  I = div_op(dx, dy)
%           I = div_op(dx, dy, wx, wy)
%
%   Input parameters:
%         dx    : Gradient along x
%         dy    : Gradient along y
%         wx    : Weights along x
%         wy    : Weights along y
%
%   Output parameters:
%         I     : Output divergence image 
%
%   Compute the 2-dimensional divergence of an image. If a cube is given,
%   it will compute the divergence of all images in the cube.
%
%   Warning: computes the divergence operator defined as minus the adjoint
%   of the gradient 
%
%   ..      div  = - grad'
%
%   .. math:: \text{div} = - \nabla^*
%
%   See also: gradient_op div_op3d div_op1d laplacian_op prox_tv

% Author: Nathanael Perraudin
% Date:   1 February 2014

if nargin > 2
    dx = dx .* conj(wx);
    dy = dy .* conj(wy);
end

I = [dx(1, :,:) ; ...
    dx(2:end-1, :,:)-dx(1:end-2, :,:) ;...
    -dx(end-1, :,:)];
I = I + [dy(:, 1,:) ,...
    dy(:, 2:end-1,:)-dy(:, 1:end-2,:) ,...
    -dy(:, end-1,:)];

end
