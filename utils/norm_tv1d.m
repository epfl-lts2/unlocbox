function y = norm_tv1d(I,w)
%NORM_TV1D 1 Dimentional TV norm
%   Usage:  y = norm_tv1d(x)
%           y = norm_tv1d(x,w)
%
%   Input parameters:
%         I     : Input data 
%         w    : Weights
%   Output parameters:
%         y     : Norm
%
%   Compute the 1-dimentional TV norm of I. If the input I is a matrix.
%   This function will compute the norm of all line and return a vector of
%   norms.
%
%   See also: norm_tv norm_tv3d


% Author: Nathanael Perraudin
% Date:   1 February 2014

if nargin>1
    dx = gradient_op1d(I,w);
else
    dx = gradient_op(I);
end    

%y = sum(temp(:));
y = sum(abs(dx),1)';

end
