function [dx ] = gradient_op1d(I, wx)
%GRADIENT_OP1D 1 Dimentional gradient operator
%   Usage:  dx = gradient_op1d(I)
%           dx = gradient_op1d(I, wx)
%
%   Input parameters:
%         I     : Input data 
%         wx    : Weights along x
%
%   Output parameters:
%         dx    : Gradient along x
%
%   Compute the 1-dimentional gradient of I. If the input I is a matrix.
%   This function will compute the gradient of all vectors and return a
%   matrix. 
%
%   See also: gradient_op gradient_op3d div_op laplacianx_op

% Author: Nathanael Perraudin
% Date:   1 February 2014

dx = [I(2:end, :)-I(1:end-1, :) ; zeros(1, size(I, 2))];

if nargin>1
    dx = dx .* wx;
end

end
