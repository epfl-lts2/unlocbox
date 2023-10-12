function reset_seed(n)
%GSP_RESET_SEED Reset the seed of the random number generator
%   Usage:  reset_seed(n);
%
%   Input parameters:
%       n   : seed
%   Ouptut parameters:
%       none
%
%   This function resets the seed

% Authors: Nathanael Perraudin, Vassilis Kalofolias
% Date  : 21 May 2014
% 


if nargin<1
    n = 0;
end


rng(n,'twister');

end

