function [weights]=test_weights(weights,warning)
%TEST_GAMMA test if the weights are corrects
%   Usage:  weights=test_weights(weights)
%           test_weights(weights)
%           weights=test_weights(weights,warning)
%           test_weights(weights,warning)
%
%   Input parameters:
%         weights : vector
%         warning: boolean
%   Output parameters:
%         weights : vector
%
%   This function test is the weights are striclty positivs
%   
%   For each element of the vector weights
%   If it is negativ, this function return an error. If it is zero
%   this function add eps to this weights and display a warning. This can be used
%   to set 0 weight to some objectiv function.
%
%   The warning parameter is a flag to display warning or not. Default=1;

% Author:  Nathanael Perraudin
% Date: Nov 2012
%

if nargin<2
   warning=1; 
end

if sum(weights<0)
    error('gamma can not be negativ!');
elseif ~logical(sum(weights(:)~=0)) && warning
    weights=weights+eps;
    fprintf(' WARNING!!! weights is 0. We add eps to weights to keep going...\n');
% else
%    weights=weights; 
end

end
    
