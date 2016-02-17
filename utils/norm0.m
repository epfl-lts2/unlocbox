function n0 = norm0(x,tol)
%NORM0 Compute the 0 norm of a vector
%   Usage: n0 = norm0(x);
%          n0 = norm0(x,tol);
%
%   Input parameters
%       x   : Vector or matrix
%       tol : Tolerance (default 1e-10)
%
%   Ouput parameters:
%       n0  : Zero norm
%
%   This function compute the zero norm of a vector or a matrix. More
%   explicitely, it count the number or non zero entries.
%

% Author: Nathanael Perraudin
% Date  : 1 July 2014


if nargin<2
    tol = 1e-10;
end

n0 = sum(abs(x(:))>tol);

end