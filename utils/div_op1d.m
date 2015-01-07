function I = div_op1d(dx, wx)
%DIV_OP1D Divergence operator in 1 dimention
%   Usage:  I = div_op1d(dx)
%           I = div_op1d(dx, wx)
%
%   Input parameters:
%         dx    : Gradient along x
%         wx    : Weights along x
%
%   Output parameters:
%         I     : Output divergence vector 
%
%   Compute the 1-dimentional divergence of a vector. If a matrix is given,
%   it will compute the divergence of all vectors in the matrix.
%
%   Warning this function compute the divergence operator defined as minus
%   the adjoint of the gradient
%
%   ..      div  = - grad'
%
%   .. math:: \text{div} = - \nabla^*
%
%   See also: gradient_op div_op3d div_op1d laplacian_op

% Author: Nathanael Perraudin
% Date:   1 February 2014

if nargin > 2
    dx = dx .* conj(wx);
end

I = [dx(1, :) ; dx(2:end-1, :)-dx(1:end-2, :) ; -dx(end-1, :)];

end
