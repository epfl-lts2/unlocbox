function set_seed(my_seed)
%SET_SEED sets the seed of the default random random generator
%   Usage:  set_seed(my_seed)
%           set_seed()
%
%   Input parameters:
%         my_seed  :  new_seed
%   
%   Set the seed of the default random random generator

% code author: Vassilis Kalofolias
% date: August 2013


%% Fix the random stream for debugging reasons
if nargin < 1
    my_seed = 0;
end

if verLessThan('matlab', '7.12.0')  % release 2011a has "rng"
    rand('twister', my_seed);
else
    rng(my_seed, 'twister');
    %rng('default');
end

end



