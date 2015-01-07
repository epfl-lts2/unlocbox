function [U, Sigma, V, r, i] = svst(A, lambda,param)

% Singular Value Soft-Thresholding (with parameter lambda)
% This is the prox of the nuclear norm i.e.,
% argmin_x  1/2 |A - x|_F^2 + lambda |x|_*
% Author: Mohammad Golbabaee
% e-mail: mohammad.golbabaei@epfl.ch
% 08.06.2011, Lausanne, Switzerland

if ~isfield(param, 's_initial'), param.s_initial = min(size(A)); end
if ~isfield(param, 'verbose'), param.verbose = 0; end

[n1,n2] = size(A);
if n1*n2 <= 100*100, SMALLSCALE = true; else SMALLSCALE = false; end

% What the best way to multiply a sparse matrix?
% Here we use Stephen Beck's code thats using fast sparse matrix
% multiplication (smvp).

% Default type.
forwardType = 1;
transposeType = 1;

% A = sparse(A);
% [forwardType, transposeType] = findBestMultiply(A,.2);
% A = full(A);

incre = 5;
s = param.s_initial;
i = 0;

if SMALLSCALE
    [U,Sigma,V] = svd(A,'econ');
else
    % Make routines for multiplying by a sparse matrix    
    At = A';
    switch forwardType
        case 1, Aforward = @(x) A*x;
        case 2, Aforward = @(x) At'*x;
        case 3, A = sparse(A); Aforward = @(x) smvp(A,x);
    end
    switch transposeType
        case 1, Atranspose = @(x) At*x;
        case 2, Atranspose = @(x) A'*x;
        case 3, At = sparse(At); Atranspose = @(x) smvp(At,x);
    end
    OK = 0;
    while ~OK
        opts = [];
        %if ~isreal(b), opts.eta = 1e-16; end
        [U,Sigma,V] = lansvd(Aforward,Atranspose,n1,n2,s,'L',opts);
        %[U,Sigma,V] = lansvd(Y,s,'L');
        OK = (Sigma(s,s) <= lambda) || ( s == min(n1,n2) );
        s = min(s + incre, min(n1,n2));
        i = i+1;
    end
end

sigma = diag(Sigma); r = sum(sigma > lambda);
U = U(:,1:r); V = V(:,1:r); sigma = sigma(1:r) - lambda; Sigma = diag(sigma);


if param.verbose==1,  fprintf('Nb_iterations %4d, rank is %2d, \n',i,r); end
