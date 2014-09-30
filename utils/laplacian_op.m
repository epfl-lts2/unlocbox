function [ L ] = laplacian_op( I )
%LAPLACIAN_OP 2 dimentional Laplacian
%   Usage:  [I] = laplacian_op( I );
%
%   Input parameters:
%         I     : Input image 
%
%   Output parameters:
%         I     : Laplacian
%
%   Compute the sum of the laplacian along x and y. This operator is
%   self-adjoint.
%
%   ..      L = I_xx + I_yy
%
%   .. math::  \mathcal{L} = I_{xx} + I_{yy}
%
%   See also: laplacianx_op laplaciany_op div_op gradient_op
%

% Author: Nathnaael Perraudin
% Date  : 13 September 2013


L=laplacianx_op(I)+laplaciany_op(I);

end

