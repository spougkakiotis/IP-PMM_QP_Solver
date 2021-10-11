function [R_NE_PS] = Reg_NE_preconditioner(A,G,A_tr,delta,R,C)
% ======================================================================================================== %
% Build Block diagonal preconditioner for the normal equations BB^T + delta I.
% -------------------------------------------------------------------------------------------------------- %
% [L_11,perm_11,L_22,perm_22] = Reg_NE_preconditioner(B,delta,k,Perm_r,Perm_c)
% 
% Input: 
%   A -> an mxn sparse matrix, (m <= n)
%   G -> an nx1 positive diagonal vector representing diag(G),
%   A_tr -> the transpose of A used to create the normal equations AGA_tr 
%   delta -> the regularization parameter,
%   R -> a logical vector indicating with true which rows are NOT dropped.
%   C -> a logical vector indicating with true which columns are NOT dropped.
%
% Output: R_NE_PS which is a struct containing all information need to build the preconditioner. 
%   Specifically:
%       .L_11 -> The lower Cholesky factor of the (1,1) block of the regularized normal equations
%       .perm_11 -> the permutation for the aforementioned Cholesky factor
%       .L_22 -> The lower Cholesky factor of an approximation of the (2,2) block of the regularized
%               normal equations
%       .perm_22 -> the permutation for the aforementioned Cholesky factor.
%       .Perm_r -> the row permutation needed for applying the preconditioner
%       .instability_11 -> if not zero, the (1,1) block is unstable.
%       .instability_22 -> if not zero, the (2,2) block is unstable.
% ________________________________________________________________________________________________________ %
    R_NE_PS = struct();
    % ==================================================================================================== %
    % Input control.
    % ---------------------------------------------------------------------------------------------------- %
    if (nargin < 4)
        fprintf("Not enough input arguments.\n");
        return;
    end
    [m, n] = size(A);  
    R_NE_PS.m = m;
    if (nargin < 5 || isempty(R))
        R = true(m,1);
    end
    if (nargin < 6 || isempty(C))
        C = true(n,1);
    end

    not_R = ~R;
    %A_1_tr = A_tr(:,not_R);
    A_1_tr = A_tr(C,not_R);
    A_1 = A_1_tr';
    A_2 = A(R,C);
    Perm_r = [find(not_R); find(R)];
    % ____________________________________________________________________________________________________ %
 
    % ==================================================================================================== %
    % Form the approximate normal equations matrix and compute the Cholesky for both blocks
    % ---------------------------------------------------------------------------------------------------- %
    R_NE_PS.k_r = nnz(not_R);
    R_NE_PS.k_c = nnz(C);
    if (R_NE_PS.k_r)
         P_11 = A_1*spdiags(G(C),0,R_NE_PS.k_c,R_NE_PS.k_c)*A_1_tr + delta.*speye(R_NE_PS.k_r);
      %  P_11 = A_1*spdiags(G,0,n,n)*A_1_tr + delta.*speye(R_NE_PS.k_r);
        [R_NE_PS.L_11,R_NE_PS.instability_11,R_NE_PS.perm_11] = chol(P_11,'lower','vector'); 
        maxpiv_11 = max(spdiags(P_11,0));
    else
        R_NE_PS.L_11 = [];
        R_NE_PS.perm_11 = [];
        R_NE_PS.instability_11 = 0;
        maxpiv_11 = 0;
    end

    P_22 = A_2*spdiags(G(C),0,R_NE_PS.k_c,R_NE_PS.k_c)*A_2' + delta.*speye(m-R_NE_PS.k_r);
    maxpiv_22 = max(spdiags(P_22,0));
    if ((~isnan(maxpiv_11) && ~isinf(maxpiv_11)) && (~isnan(maxpiv_11) && ~isinf(maxpiv_11)))
        R_NE_PS.maxpiv_11 = maxpiv_11;   R_NE_PS.maxpiv_22 = maxpiv_22;
        if (nnz(C))
            [R_NE_PS.L_22,R_NE_PS.instability_22,R_NE_PS.perm_22] = chol(P_22,'lower','vector');
            if (R_NE_PS.instability_11 ~= 0 || R_NE_PS.instability_22 ~= 0)
                fprintf("The matrix is badly conditioned. Abort and recompute.\n");
            end  
        else
            R_NE_PS.L_22 = speye(m-R_NE_PS.k_r);
            R_NE_PS.instability_22 = 0;
            R_NE_PS.perm_22 = 1:m;
        end
    end
    if (R_NE_PS.instability_11 ~= 0 || R_NE_PS.instability_22 ~= 0)
        fprintf("The matrix is badly conditioned. Abort and recompute.\n");
    end
    R_NE_PS.Perm_r = sparse((1:m)',Perm_r,ones(m,1));
    % ____________________________________________________________________________________________________ %  
end


