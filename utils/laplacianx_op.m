function [ Lx ] = laplacianx_op( I )
%LAPLACIANX_OP dimentional Laplacian
%   Usage:  [Lx] = laplacianx_op( I );
%
%   Input parameters:
%         I     : Input image 
%
%   Output parameters:
%         Lx    : Laplacian along x
%
%   Compute the sum of the laplacian along x. This operator is 
%   self-adjoint.
%
%   ..      Lx = I_xx
%
%   .. math::  \mathcal{L}_x = I_{xx}
%
%   See also: laplacian_op laplaciany_op div_op gradient_op
%

% Author: Nathnaael Perraudin
% Date  : 13 September 2013


dx = [I(2:end, :,:)-I(1:end-1, :,:) ; zeros(1, size(I, 2),size(I, 3))];


Lx = [dx(1, :,:) ; dx(2:end-1, :,:)-dx(1:end-2, :,:) ; -dx(end-1, :,:)];


end
