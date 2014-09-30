function ninf1 = norm_linf1(x, g_d,g_t, winf,w1)
%NORM_Linf1 Linf1 mixed norm
%   Usage:  ninf1 = norm_linf1(x);
%           ninf1 = norm_linf1(x, g_d,g_t);
%           ninf1 = norm_linf1(x, g_d,g_t, winf,w1);
%
%   Input parameters:
%         x     : Input data 
%         g_d   : group vector 1
%         g_t   : group vector 2
%         winf  : weights for the sup norm (default 1)
%         w1    : weights for the one norm (default 1)
%   Output parameters:
%         y     : Norm
%
%   `norm_linf1(x, g_d,g_t, w2,w1)` returns the norm Linf1 of x. If x is a
%   matrix the sup norm will be computed over the lines (2nd dimention) and
%   the one norm will be computed over the rows (1st dimention). In this
%   case, all other argument are not necessary.
%
%   ..  ninf1 = || x ||_inf1 = sum_j (max_i |x(i,j)|)
%
%   .. math::  \| x \|_{\infty 1} = \sum_j \left| \max_i |x(i,j)|  \right| 
%
%   'norm_linf1(x)' with x a row vector is equivalent to norm(x,1) and
%   'norm_linf1(x)' with x a line vector is equivalent to max(abs(x))
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
%   See also: norm_l21 norm_tv


% Author: Nathanael Perraudin
% Date: October 2011
% Testing: test_mixed_sparsity


% Optional input argumentsX
if nargin<2, g_d = 1:numel(x); end
if nargin<3, 
    if numel(x) == size(x,1)*size(x,2); % matrix case
        g_t = size(x,2)*ones(1,size(x,1)); 
    else
        g_t = ones(1,numel(x)); 
    end
end
if nargin<5, w1=ones(size(g_t,2),size(g_t,1)); end
if nargin<4, winf=ones(numel(x),size(g_t,1)); end



% overlapping groups
if size(g_d,1)>1
    n21=0;
    for ii=1:size(g_d,1);
        n21 = n21 + norm_linf1(x,g_d(ii,:),g_t(ii,:), ...
            winf(:,ii),w1(:,ii));
    end
else % non overlapping groups

    l=length(g_t);

    % Compute the norm
    if max(g_t)==min(g_t),
        X=transpose(x);
        X=X(g_d);
        X=transpose(reshape(X,numel(x)/l,l));
        Winf=transpose(reshape(winf(g_d),numel(x)/l,l));
        normsup=max(abs(Winf.*X),[],2);
        ninf1=sum(w1.*normsup);
    else % group of different size

        ninf1=0;
        indice=0;
        X=x;
        X=X(g_d);
        Winf=winf;
        Winf=Winf(g_d);
        for ii=1:l

            ninf1=ninf1+w1(ii)*norm(Winf(indice+1:indice+g_t(ii)).*X(indice+1:indice+g_t(ii)),Inf);

            indice=indice+g_t(ii);
        end

    end

end


end
