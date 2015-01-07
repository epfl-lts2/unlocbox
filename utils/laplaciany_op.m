function [ Ly ] = laplaciany_op( I )
%LAPLACIANY_OP dimentional Laplacian
%   Usage:  [Ly] = laplaciany_op( I );
%
%   Input parameters:
%         I     : Input image 
%
%   Output parameters:
%         Ly    : Laplacian along y
%
%   Compute the sum of the laplacian along y. This operator is
%   self-adjoint.
%
%   ..      Ly = I_yy
%
%   .. math::  \mathcal{L}_y = I_{yy}
%
%   See also: laplacian_op laplacianx_op div_op gradient_op
%

% Author: Nathnaael Perraudin
% Date  : 13 September 2013


dy = [I(:, 2:end,:)-I(:, 1:end-1,:) , zeros(size(I, 1), 1,size(I, 3))];

Ly =  [dy(:, 1,:) , dy(:, 2:end-1,:)-dy(:, 1:end-2,:) , -dy(:, end-1,:)];

end
