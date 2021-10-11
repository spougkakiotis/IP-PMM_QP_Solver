function [R] = find_k_drop_rows(k, A, m, mode)
% ================================================================================ %
% [R] = find_k_drop_rows(k, row, m)
% -------------------------------------------------------------------------------- %
% Input:
%   k: number of most dense rows to be tracked.
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
    if (nargin < 4 || isempty(mode))
        mode = "Inf_norm";
    end
    if (mode == "Density")
        [row,~,~] = find(A);
        row_nnzs = zeros(m,1);
        for i = 1:m
            row_nnzs(i) = nnz(find(row == i));
        end
        [elems,dense_rows] = maxk(row_nnzs,k);   % Indexes of the k_r most dense rows.
        elems
    else
        sum_rows = sum(A,2);
        [~,dense_rows] = mink(sum_rows,k);   % Indexes of the k_r most dense rows.
    end
    R(dense_rows) = false;    
end

