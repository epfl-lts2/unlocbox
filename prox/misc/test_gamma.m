function [ stop ]=test_gamma(gamma)
%TEST_GAMMA test if gamma is correct
%   Usage:  stop = test_gamma(gamma)
%           test_gamma(gamma)
%
%   Input parameters:
%         gamma : number
%   Output parameters:
%         stop  : boolean
%
%   This function test is gamma is stricly positive
%   
%   If gamma is negativ, this function return an error. If gamma is zero
%   this function, set stop to 1.
%
%

% Author:  Nathanael Perraudin
% Date: February 2012
%


if gamma<0
    error('gamma can not be negativ!');
% elseif (gamma==0) && warning
%     gamma=gamma+eps;
%     fprintf(' WARNING!!! gamma is 0. We add eps to gamma to keep going...\n');
% else
%   % gamma = gamma; 
end
    
if gamma==0
    stop = 1;
else
    stop = 0;
end

end
