function [U,S,V] = svdsecon(X,k)
%SVDSECON Fast svds when n<<m
%   Usage: [U,S,V] = svdsecon(X,k);
%
%   Input parameters:
%         X     : Input data (n x m)
%         k     : Number of singular values
%   Output parameters:
%         U     : Left singular vectors
%         S     : Singular values
%         U     : Right signular vectors
%
%   This function is an acceleration of svds. It is particularly efficient
%   when n<<m
%
%   See also: svdecon

    [m,n] = size(X);

    if  m <= n
        X2 = X*X';
        [U,E] = eigs(X2,k);

        [e,ix] = sort(abs(diag(E)),'descend');
        U = U(:,ix);    

        V = X'*U;
        s = sqrt(e);
        
        V = bsxfun(@times, V, 1./s');
        S = diag(s);
    else
        X2 = X'*X; 
        [V,E] = eigs(X2,k);

        [e,ix] = sort(abs(diag(E)),'descend');
        V = V(:,ix);    

        U = X*V; 
        s = sqrt(e);
        U = bsxfun(@times, U, 1./s');
        S = diag(s);
    end
end