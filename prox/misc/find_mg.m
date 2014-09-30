function [ Kw, ny, Mg, zo,wo] = find_mg( z,lambda,w )
%FIND_MG find the MG other stuff for the prox_l12
%
%   In the litterature, they are mistakes in the formulas. As consquence,
%   this function is a mixed of the formulas found in those 3 papers. In the
%   first paper, the formulas seems correct. The easier to read is most
%   likely the third one.
%
%   Cheers!
%
%   References: kowalski2009sparse kowalski2013social kowalski2009sparsity


if nargin<3
    w = ones(size(z));
end


r = abs(z./w);

[Nelg, Ng] = size(r);

[r, ind] = sort(r,1,'descend');

zo = zeros(Nelg,Ng);
wo = zeros(Nelg,Ng);

for jj = 1 : Ng
    zo(:,jj) = z(ind(:,jj),jj);
    wo(:,jj) = w(ind(:,jj),jj);
end

Mg = zeros(Ng,1);
Kw = zeros(Ng,1);
ny = zeros(Ng,1);

for jj = 1 : Ng

    for ii=1:Nelg-1
        temp = lambda * ...
            sum( wo(1:(ii+1),jj).^2 .* ( r(1:(ii+1),jj) - r((ii+1),jj) ) )...
            - r((ii+1),jj);
        if temp >= 0
            Mg(jj) = ii;
            Kw(jj) = sum( wo(1:ii,jj).^2 );
            ny(jj) = norm( zo(1:ii,jj),1);
            break;
        end
    end
    
    % handle exception
    if Mg(jj)==0
        Mg(jj) = 1;
        Kw(jj) = wo(1,jj).^2;
        ny(jj) = abs(zo(1,jj));
    end

end




end
