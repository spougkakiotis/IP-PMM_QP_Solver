function [R_NE_PS] = Reg_NE_preconditioner(B,delta,k,Perm_r,Perm_c)
% ======================================================================================================== %
% Build Block diagonal preconditioner for the normal equations BB^T + delta I.
% -------------------------------------------------------------------------------------------------------- %
% [L_11,perm_11,L_22,perm_22] = Reg_NE_preconditioner(B,delta,k,Perm_r,Perm_c)
% 
% Input: 
%   B -> an mxn sparse matrix, (m <= n)
%   delta -> the regularization parameter,
%   k -> i) either an integer k in [1,m] indicating the number of rows and columns dropping 
%        ii) or a vector of size 2x1, k = [k_r,k_c] in k_r in [1,m], k_c in [1,n] indicating the 
%            repsective number of rows and columns to be dropped.
%   Perm_r -> a row permutation vector, whose first k_r (or k) elements point to the k_r most dense 
%             rows of B. (can be empty, or missing if Perm_c is also missing)
%   Perm_c -> a column permutation vector, whose first k_c (or k) elements point to the k_c most dense
%             columns of B. (can be empty, or missing)
%
% Output: R_NE_PS which is a struct containing all information need to build the preconditioner. 
%   Specifically:
%       .L_11 -> The lower Cholesky factor of the (1,1) block of the regularized normal equations
%       .perm_11 -> the permutation for the aforementioned Cholesky factor
%       .L_22 -> The lower Cholesky factor of an approximation of the (2,2) block of the regularized
%               normal equations
%       .perm_22 -> the permutation for the aforementioned Cholesky factor.
%       .Perm_r -> the row permutation needed for applying the preconditioner
%       .Perm_c -> the column permutation needed for applying the preconditioner
%       .instability_11 -> if not zero, the (1,1) block is unstable.
%       .instability_22 -> if not zero, the (2,2) block is unstable.
% ________________________________________________________________________________________________________ %
    R_NE_PS = struct();
    % ==================================================================================================== %
    % Input control.
    % ---------------------------------------------------------------------------------------------------- %
    if (nargin < 3)
        fprintf("Not enough input arguments.\n");
        return;
    end
    [m,n] = size(B);
    R_NE_PS.m = m;
    R_NE_PS.n = n;
    compute_perm_r = false; compute_perm_c = false;
    L_11 = [];  perm_11 = [];   L_22 = [];  perm_22 = [];
    if (size(k,1) == 1)
        k_r = k;    k_c = k;
    elseif (size(k,1) == 2)
        k_r = k(1); k_c = k(2);
    else
        fprintf("Incorrect input argument for number of dropping rows or columns.\n");
        return;
    end
        
    if (nargin < 4 || isempty(Perm_r))
        compute_perm_r = true;
        Perm_r = 1:m;
    end
    if (nargin < 5 || isempty(Perm_c))
        compute_perm_c = true;
        Perm_c = 1:n;
    end
    % ____________________________________________________________________________________________________ %
 
    % ==================================================================================================== %
    % If not provided, compute row and column permutations that bring the k_r most dense row and 
    % k_c most dense columns of B at the top.
    % ---------------------------------------------------------------------------------------------------- %
    [row,col,v] = find(B);
    nnz_B = size(row,1);
    if (compute_perm_r || compute_perm_c)
        dense_rows = [];    dense_cols = [];
    end
    if (compute_perm_r)
        row_nnzs = zeros(m,1);
        for i = 1:m
            row_nnzs(i) = nnz(find(row == i));
        end
        [~,dense_rows] = maxk(row_nnzs,k_r);   % Indexes of the k_r most dense rows.
        for j = 1:k_r
            tmp = Perm_r(j);
            Perm_r(j) = dense_rows(j);
            Perm_r(dense_rows(j)) = tmp;              % The row permutation matrix
        end
    end
    if (compute_perm_c)
        col_nnzs = zeros(m,1);
        for i = 1:n
            col_nnzs(i) = nnz(find(col == i));
        end
        [~,dense_cols] = maxk(col_nnzs,k_c);   % Indexes of the k_r most dense rows. 
        for j  = 1:k_c
            tmp = Perm_c(j);
            Perm_c(j) = dense_cols(j);
            Perm_c(dense_cols(j)) = tmp;              % The column permutation matrix
        end
    end
    % ____________________________________________________________________________________________________ %

    % ==================================================================================================== %
    % Use the perumation vector to permute the matrix in the sparse format.
    % ---------------------------------------------------------------------------------------------------- %
    row_permuted = zeros(nnz_B,1);
    col_permuted = zeros(nnz_B,1);
    for k = 1:nnz_B
        row_permuted(k) = Perm_r(row(k));             % Permuting rows in sparse format (i,j,v) idxs
        col_permuted(k) = find(Perm_c == col(k));     % Permuting columns in sparse format.
    end
    % ____________________________________________________________________________________________________ %

    % ==================================================================================================== %
    % Using the permuted sparse representation, create each of the three blocks of interest directly
    % from their sparse indexes. Note that the (2,1) block will not be used in the approximation.
    % ---------------------------------------------------------------------------------------------------- %
    col_11 = find(col_permuted <= k_c);
    row_11 = find(row_permuted <= k_r);
    rc_11 = intersect(col_11,row_11);
    B_11 = sparse(row_permuted(rc_11),col_permuted(rc_11),v(rc_11),k_r,k_c); % Create (1,1) block
    col_12 = find(col_permuted > k_c);
    rc_12 = intersect(col_12,row_11);
    B_12 = sparse(row_permuted(rc_12),col_permuted(rc_12)-k_c,v(rc_12),k_r,n-k_c);  % Create (1,2) block
    col_22 = find(col_permuted > k_c);
    row_22 = find(row_permuted > k_r);
    rc_22 = intersect(col_22,row_22);
    B_22 = sparse(row_permuted(rc_22)-k_r,col_permuted(rc_22)-k_c,v(rc_22),m-k_r,n-k_c); % (2,2) block.
    % ____________________________________________________________________________________________________ %
    
    % ==================================================================================================== %
    % Form the approximate normal equations matrix and compute the Cholesky for both blocks
    % ---------------------------------------------------------------------------------------------------- %
    M_11 = B_11*B_11' + B_12*B_12' + spdiags(delta.*ones(k_r,1),0,k_r,k_r);
    M_22 = B_22*B_22' + spdiags(delta.*ones(m-k_r),0,m-k_r, m-k_r);
    [R_NE_PS.L_11,R_NE_PS.instability_11,R_NE_PS.perm_11] = chol(M_11,'lower','vector'); 
    [R_NE_PS.L_22,R_NE_PS.instability_22,R_NE_PS.perm_22] = chol(M_22,'lower','vector');
    if (R_NE_PS.instability_11 ~= 0 || R_NE_PS.instability_22 ~= 0)
        fprintf("The matrix is badly conditioned. Abort and recompute.\n");
    end
    R_NE_PS.Perm_c = sparse((1:n)',Perm_c',ones(n,1));
    R_NE_PS.Perm_r = sparse((1:m)',Perm_r',ones(m,1));
    R_NE_PS.k_r = k_r;  R_NE_PS.k_c = k_c;
    % ____________________________________________________________________________________________________ %  
end


