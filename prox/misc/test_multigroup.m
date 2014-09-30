function []=test_multigroup(x,g_d,g_t)
%TEST_MULTIGROUP test if the parameter g_d and g_t are correct
%   Usage:  test_multigroup(x,g_d,g_t)
%
%   Input parameters:
%         x       : vector
%         g_d     : numerics
%         g_t     : numerics
%   Output parameters:
%
%   This function 
%

% Author:  Nathanael Perraudin
% Date: Nov 2012
%


L=numel(x); % lenght of the signal

if size(g_d,2)~=L
   fprintf(' WARNING!!! The number of collum of gd should be the size of the signal and it is not.\n'); 
end

if sum(g_t,2)~=size(g_d,2)
   error( ' The sum of the group size should be equal to the total number of element grouped. Check if g_d,g_t are row vectors!');
end


