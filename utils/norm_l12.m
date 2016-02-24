function n12 = norm_l12(x, g_d,g_t,w1,w2)
%NORM_L12 L12 mixed norm
%   Usage:  n12 = norm_l12(x);
%           n12 = norm_l12(x, g_d,g_t);
%           n12 = norm_l12(x, g_d,g_t, w1,w2);
%
%   Input parameters:
%         x     : Input data 
%         g_d   : group vector 1
%         g_t   : group vector 2
%         w1    : weights for the one norm (default 1)
%         w2    : weights for the two norm (default 1)
%   Output parameters:
%         y     : Norm
%
%   `norm_l12(x, g_d, g_t, w1, w2)` returns the norm L21 of x. If $x$ is a
%   matrix, the norm will be computed as follow: 
%
%   ..  n12 = || x ||_12 = ( sum_j | sum_i |x(i,j)| |^2 )^(1/2) 
%
%   .. math::  \| x \|_{1,2}^2 = \sqrt{ \sum_j \left| \sum_i |x(i,j)| \right|^2 }
%
%   In this case, all other argument are not necessary.
%
%   'norm_l12(x)' with x a row vector is equivalent to norm(x,2) and
%   'norm_l12(x)' with x a line vector is equivalent to norm(x,1)
%
%   For fancy group, please provide the groups vectors.
%
%   `g_d`, `g_t` are the group vectors. `g_d` contain the indices of the
%   element to be group and `g_t` the size of different groups.
%       
%   Example: 
%                x=[x1 x2 x3 x4 x5 x6] 
%                Group 1: [x1 x2 x4 x5] 
%                Group 2: [x3 x6]
%
%   Leads to 
%           
%               => g_d=[1 2 4 5 3 6] and g_t=[4 2]
%               Or this is also possible
%               => g_d=[4 5 3 6 1 2] and g_t=[2 4]   
%
%   This function works also for overlapping groups.
%
%   See also: norm_linf1 norm_l21
%
%   References: kowalski2013social



% Author: Nathanael Perraudin
% Date: Mai 2014
% Testing: test_mixed_sparsity


% Optional input arguments
if nargin<2, g_d = 1:numel(x); end
if nargin<3, 
    if numel(x) == size(x,1)*size(x,2); % matrix case
        g_t = size(x,2)*ones(1,size(x,1)); 
    else
        g_t = ones(1,numel(x)); 
    end
end
if nargin<5, w2=ones(size(g_t,2),size(g_t,1)); end
if nargin<4, w1=ones(numel(x),size(g_t,1)); end


% overlapping groups
if size(g_d,1)>1
    n12=0;
    for ii=1:size(g_d,1);
        n12 = n12 + norm_l12(x,g_d(ii,:),g_t(ii,:), ...
            w1(:,ii),w2(:,ii));
    end
else % non overlapping groups


    l=length(g_t);

    % Compute the norm
    if max(g_t)==min(g_t) % group of the same size
        X = transpose(x);
        X = X(g_d);
        X = transpose(reshape(X,numel(x)/l,l));

        W1 = transpose(reshape(w1(g_d),numel(x)/l,l));
        normX1 = sum(abs(W1.*X),2);
        n12 = norm(w2.*normX1,2);
    else % group of different size

        n12 = zeros(l,1);
        indice = 0;
        X = x(:);
        X = X(g_d);
        W1 = w1;
        W1 = W1(g_d);
        for ii=1:l

            n12(ii) =  w2(ii) * norm(W1(indice+1:indice+g_t(ii)) ...
                         .*X(indice+1:indice+g_t(ii)),1);
            indice = indice+g_t(ii);
        end

        n12 = norm(n12,2);
    end
    
end

end
