function n = back_perm(k)
%BACK_PERM(k)
%   Usage: n = back_perm(k)
%
%   Input parameters:
%         k     : row vectors of permutation (Can be a matrix)
%   Output parameters:
%         n     : Permuted vector
%
%   Warning! k has to be a row vector.

% Author:  Nathanael Perraudin
% Date: October 2011
%
n=zeros(size(k));

for j=1:size(k,2)
for p=1:size(k,1)
    n(k(p),j)=p;
end
end

end
