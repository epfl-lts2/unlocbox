function w = squareform_sp(w)
%SQUAREFORM_SP Sparse counterpart of matlab's squareform
%   Usage: w = squareform_sp(W);
%
%   Input parameters:
%       w: sparse vector with n(n-1)/2 elements OR
%       W: matrix with size [n, n] and zero diagonal
%
%   Output parameters:
%       W: matrix form of input vector w OR
%       w: vector form of input matrix W
%
%   This function is to be used instead of squareform.m when the matrix W
%   or the vector w is sparse. For large scale computations, e.g. for
%   learning the graph structure of a big graph it is necessary to take
%   into account the sparsity.
%
%   Example::
%
%       B = sprand(8, 8, 0.1);
%       B = B+B';
%       B(1:9:end) = 0;
%       b = squareform_sp(B);
%       Bs = squareform_sp(b);
%
%
%   See also: sum_squareform zero_diag

%   Date: December 2015
%   Author: Vassilis Kalofolias
%   Testing: test_squareform

% if input is not sparse it doesn't make sense to use this function!
if not(issparse(w)) && nnz(w)/numel(w) > 1/10
    w = squareform(w);
%     warning('Use standard squareform for non sparse vector! Full version returned.');
    return
end


if isvector(w)
    %% VECTOR -> MATRIX
    l = length(w);
    n = round((1 + sqrt(1+8*l))/2);
    % check input
    if not(l == n*(n-1)/2)
        error('Bad vector size!');
    end
    
    %% TODO: make more efficient! (Maybe impossible)
    if iscolumn(w)
        [ind_vec, ~, s] = find(w);
    else
        [~, ind_vec, s] = find(w);
    end
    
    num_nz = length(ind_vec);
    
    % indices inside the matrix
    ind_i = zeros(num_nz, 1);
    ind_j = zeros(num_nz, 1);
    
    curr_row = 1;
    offset = 0;
    % length of current row of matrix, counting from after the diagonal
    len = n - 1;
    for ii = 1 : length(ind_vec)
        ind_vec_i = ind_vec(ii);
        % if we change row, the vector index is bigger by at least the
        % length of the line + the offset.
        while(ind_vec_i > (len + offset))
            offset = offset + len;
            len = len - 1;
            curr_row = curr_row + 1;
        end
        ind_i(ii) = curr_row;
        ind_j(ii) = ind_vec_i - offset + (n-len);
    end
    
    % for the lower triangular part just add the transposed matrix
    %w = sparse(ind_i, ind_j, s, n, n) + sparse(ind_j, ind_i, s, n, n);
    w = sparse([ind_i; ind_j], [ind_j; ind_i], [s(:); s(:)], n, n);
    
else
    %% MATRIX -> VECTOR
    % first checks
    [m, n] = size(w);
    if m ~= n || ~all(diag(w)==0)
        error('Matrix has to be square with zero diagonal!');
    end
    
    [ind_i, ind_j, s] = find(w);
    % keep only upper triangular part
    ind_upper = ind_i < ind_j;
    ind_i = ind_i(ind_upper);
    ind_j = ind_j(ind_upper);
    s = s(ind_upper);
    % compute new (vector) index from (i,j) (matrix) indices
    new_ind = ind_j + (ind_i-1)*n - ind_i.*(ind_i+1)/2;
    w = sparse(new_ind, 1, s, n*(n-1)/2, 1);
end    

end