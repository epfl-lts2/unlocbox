function [U,S,V] = svdsecon(X,k)
%SVDECON Fast svds when n<<m

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