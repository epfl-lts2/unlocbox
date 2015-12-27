function [U,S,V] = svdecon(X)
%SVDECON Fast svds when n<<m
%   Usage: [U,S,V] = svdecon(X);
%
%   Input parameters:
%         X     : Input data (n x m)
%   Output parameters:
%         U     : Left singular vectors
%         S     : Singular values
%         U     : Right signular vectors
%
%   This function is an acceleration of svd. It is particularly efficient
%   when n<<m
%
%   See also: svdsecon

    [m,n] = size(X);

    if  m <= n
        X2 = X*X';
        [U,E] = eig(X2);

        [e,ix] = sort(abs(diag(E)),'descend');
        U = U(:,ix);    

        V = X'*U;
        s = sqrt(e);
        
        V = bsxfun(@times, V, 1./s');
        S = diag(s);
    else
        X2 = X'*X; 
        [V,E] = eig(X2);

        [e,ix] = sort(abs(diag(E)),'descend');
        V = V(:,ix);    
        %% TODO: compute U = X*(V/S) and NOT U = (X*V)/S
        U = X*V; 
        s = sqrt(e);
        U = bsxfun(@times, U, 1./s');
        S = diag(s);
    end
end