function [w] = RNE_precond_operator(x,R_NE_PS)
% ===================================================================================================================== %
% Preconditioner Operator:
% --------------------------------------------------------------------------------------------------------------------- %
% [w] = RNE_precond_operator(x,R_NE_PS) takes as an input the struct containing the NE preconditioner blocks, as well as 
% a vector of size m, and returns the matrix-vector product of the inverse preconditioner by this vector.
% _____________________________________________________________________________________________________________________ % 

        w = R_NE_PS.Perm_r*x;                      % Apply the permutation P_r (row permutation of preconditioner)
        w_1 = w(1:R_NE_PS.k_r,1);
        w_2 = w(R_NE_PS.k_r + 1:R_NE_PS.m,1);
        w_1 = w_1(R_NE_PS.perm_11);                 % Cholesky permutation for the (1,1) block.
        w_2 = w_2(R_NE_PS.perm_22);
        w_1 = (R_NE_PS.L_11'\(R_NE_PS.L_11\w_1));   % Backsolve of (1,1) block 
        w_2 = (R_NE_PS.L_22'\(R_NE_PS.L_22\w_2));   % Backsolve of (2,2) block
        w_1(R_NE_PS.perm_11) = w_1;                 % Undo Cholesky permutations
        w_2(R_NE_PS.perm_22) = w_2;
        w = [w_1; w_2];
        w = R_NE_PS.Perm_r'*w;                       % Undo row permutation of preconditioner.
   
end

