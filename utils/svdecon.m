function [U,S,V] = svdecon(X)
%SVDECON Fast svd when n<<m

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