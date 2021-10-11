function [R, C] = find_row_and_col_drop(k, A, A_tr, m, n, dens_c, dens_r)
% ================================================================================ %
% [R] = find_k_drop_rows(k, row, m)
% -------------------------------------------------------------------------------- %
% Input:
%   k: maximum number of most dense columns and rows to be tracked.
%   A: is a sparse matrix
%   m: the number of rows o A.
%   mode: strategy for choosing which variables to drop.
%       "Density" -> drop the k most dense rows,
%       "Inf_norm" -> drop the k rows with the largest 1-norm
% Output:
%   R: a boolean vector of size m, where R(J) = 0, |J| = k, and J corresponds
%      to indices of the k most dense rows.
% _______________________________________________________________________________ %
    R = true(m,1);
    C = true(n,1);
    if (nargin < 6 || isempty(dens_c))
        dens_c = 0.1;
    end
    if (nargin < 7 || isempty(dens_r))
        dens_r = 0.4;
    end
    row_nnzs = zeros(m,1);
    col_nnzs = zeros(n,1);
    dense_rows = [];
    dense_cols = [];
    % =========================================================================== %
    % Find all columns that have at least 'dens' density, 
    % where    'dens' in (0,1].
    % Store their indexes in dense_cols
    % --------------------------------------------------------------------------- %
    [~,col,~] = find(A);
    total_nnzs = size(col,1);
    j = 1;
    for i = 1:n
        while ((j <= total_nnzs) && (col(j) == i))
            col_nnzs(i) = col_nnzs(i) + 1;
            j = j + 1;
        end
        if (col_nnzs(i) > dens_c*m)
            dense_cols = [dense_cols; i];
        end
    end
    k_c = size(dense_cols,1);
    if (k_c >= k) 
        [~,I] = sort(col_nnzs(dense_cols),'descend');
        dense_cols = dense_cols(I); % Sort them.
        C(dense_cols(1:k)) = false;
        fprintf("We will drop %d dense columns.\n",k);
        fprintf("We will ""drop"" 0 dense rows.\n");
        return;
    elseif (k_c > 0)
        C(dense_cols) = false;
    end
    fprintf("We will drop %d dense columns.\n",k_c);
    % __________________________________________________________________________ %
        
    % ========================================================================== %
    % Find all rows that have at least 'dens' density,
    % Store their indexes in dense_rows
    % -------------------------------------------------------------------------- %
    [~,row,~] = find(A_tr);
    j = 1;
    for i = 1:m
        while ((j <= total_nnzs) && (row(j) == i))
            row_nnzs(i) = row_nnzs(i) + 1;
            j = j + 1;
        end
        if (row_nnzs(i) > dens_r*n)
            dense_rows = [dense_rows; i];
        end
    end
    k_r = size(dense_rows,1);
    if (k_r >= k - k_c) 
        [~,I] = sort(row_nnzs(dense_rows),'descend');
        dense_rows = dense_rows(I); % Sort them.
        R(dense_rows(1:(k - k_c))) = false;
        fprintf("We will ""drop"" %d dense rows.\n",k-k_c);
        return;
    elseif (k_r > 0)
        R(dense_rows) = false;
    end
    fprintf("We will ""drop"" %d dense rows.\n",k_r);
    % _________________________________________________________________________ %
end

