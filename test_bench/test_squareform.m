function n_errors = test_squareform

n_errors = test_squareform_sp + test_sum_squareform;


end




function n_errors = test_sum_squareform

n_errors = 0;

% tolerance 
tol = 1e-10;
% size of problem to be solved
n = 10;

% create a matrix suitable for squareform
w = pdist(randn(10))';
W = squareform(w);

% create operators to be tested
[S, St] = sum_squareform(n);

% check sum
s = sum(W)';
s2 = S*w;

if norm(s - s2)/norm(s) > tol
    n_errors = n_errors + 1;
end

%% check the case including a mask
% create random mask
mask =  -1+randi(2, n*(n-1)/2, 1);

[S2, S2t] = sum_squareform(n, mask);
% if mask is given, the resulting S are the ones we would get with the
% following operations (but memory efficiently):
[ind_i, ~, ~] = find(mask(:));
% get rid of the columns of S corresponding to zeros in the mask
S3 = S(:, ind_i);
S3t = St(ind_i, :);

if norm(S2 - S3, 'fro') / norm(S3, 'fro') > tol
    n_errors = n_errors + 1;
end

if norm(S2t - S3t, 'fro') / norm(S3t, 'fro') > tol
    n_errors = n_errors + 1;
end


end


function n_errors = test_squareform_sp

n_errors = 0;

% tolerance 
tol = 1e-10;
% size of problem to be solved
n = 50;

% create a matrix suitable for squareform
w = pdist(randn(n))';
% sparsify it
w = w .* (rand(size(w)) < 0.1);

% this should also be sparse
W = squareform(w);

%% create the sparse cases
Ws = squareform_sp(sparse(w));
% transposing should not change the result
Ws2 = squareform_sp(sparse(w)');
ws = squareform_sp(sparse(Ws));


if not(issparse(ws)) || not(issparse(Ws))
    n_errors = n_errors + 1;
end

if not(all(Ws(:) == Ws2(:)))
    n_errors = n_errors + 1;
end
    

if norm(Ws - W, 'fro') / norm(W, 'fro') > tol
    n_errors = n_errors + 1;
end

if norm(ws - w, 'fro') / norm(w, 'fro') > tol
    n_errors = n_errors + 1;
end




end

