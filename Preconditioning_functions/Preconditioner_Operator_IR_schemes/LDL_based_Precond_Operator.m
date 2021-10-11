function [w] = LDL_based_Precond_Operator(x,PS)
% ===================================================================================================================== %
% Preconditioner Operator:
% --------------------------------------------------------------------------------------------------------------------- %
% [w] = LDL_based_Precond_Operator(x,PS) takes as an input the struct containing the preconditioner blocks, as well as 
% a vector of size n or n+m, and returns the matrix-vector product of the inverse preconditioner by this vector. 
% The latter can be approximated via an Iterative Refinement approach.
% Three LDL-based preconditioners are supported:
%       -> Normal equations approximation via LDL + pcg
%       -> Augmented system approximation + minres (via block diagonal preconditioning),
%       -> Augmented system approximation + minres (via LDL decomposition of the saddle point system).
% _____________________________________________________________________________________________________________________ % 
     for i = 1:PS.iter_refinement_maxit
        if ((size(x,1) > PS.m) && PS.minres_setting == "augmented_system")    % Use L|D|L^T preconditioner (NO IR possible!).
            rhs = x;
            D_M = abs(PS.D_M);
            w = Preconditioner_backsolve(PS.L_M,PS.Pinv,PS.P,rhs,D_M);
        else
            if (i == 1)
                D_M = PS.D_M;
                if (size(x,1) == PS.m)
                    rhs_m = [zeros(PS.nnz_B,1); x];
                else    % This methodology only works if we guarantee that D_M is diagonal.
                    rhs_n = x(1:PS.n,1);
                    rhs_m = [zeros(PS.nnz_B,1); x(PS.n+1:PS.n + PS.m,1)];
                end 
            else
                rhs_m(1:PS.nnz_B,1) = -PS.stability_regularizer.*w_m(1:PS.nnz_B);
                if (size(x,1) == PS.m)
                    rhs_m(PS.nnz_B+1:end,1) = x + PS.stability_regularizer.*w_m(PS.nnz_B+1:end,1);
                else
                    rhs_m(PS.nnz_B+1:end,1) = x(PS.n+1:PS.n + PS.m,1) + PS.stability_regularizer.*w_m(PS.nnz_B+1:end,1);   
                    rhs_n = x(1:PS.n,1) + PS.stability_regularizer.*w_n;
                end
            end
            if (size(x,1) ~= PS.m)
                w_n = rhs_n(PS.P_Q,1);
                w_n = PS.L_Q'\(PS.L_Q\w_n);
                w_n = w_n(PS.P_Qinv,1);
            end
            w_m = Preconditioner_backsolve(PS.L_M,PS.Pinv,PS.P,rhs_m,D_M);
            if (i == PS.iter_refinement_maxit)
                w_m = w_m(PS.nnz_B+1:end,1);
                if (size(x,1) ~= PS.m)
                    w = [w_n; w_m];
                else
                    w = w_m;
                end
            end
        end
     end
end