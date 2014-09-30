function n = norm_sumg(x, G, w)
%NORM_SUMG 2 Dimentional TV norm
%   Usage:  y = norm_sumg(x, G);
%           y = norm_sumg(x, G, w);
%
%   Input parameters:
%         x     : Input data (vector)
%         G     : The structure array of norm operator: 
%         w     : Weights (default 1)
%   Output parameters:
%         n     : Norm
%
%   `n = norm_sumg(x, G, w)` returns the sum of the norm x given in
%   the structure array G. The norm can be weighted using the parameter
%   `weights`.
%
%   See also: prox_sumg

% Author:  Nathanael Perraudin
% Date: October 2011
%


% Optional input arguments
if nargin<2, error('No input functions!'); end
if nargin<3, w=ones(length(G),1); end


% Compute the norm
n=0;

for ii=1:length(G)
    n=n+w(ii)*G{ii}.eval(x);
end

end
