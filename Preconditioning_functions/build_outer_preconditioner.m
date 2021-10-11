function [PS] = build_outer_preconditioner(A,A_tr,Q,NS,droptol,mu,R,C,double_approximation,set_proximal_inner_IR)
% ======================================================================================================== %
% Build Preconditioner with matrix E: An approximation of the Schur complement.
% -------------------------------------------------------------------------------------------------------- %
% [PS] = build_outer_preconditioner(A,A_tr,Q,NS,droptol,mu,R,C,double_approximation,set_proximal_inner_IR)
% ________________________________________________________________________________________________________ %

    PS = struct();                    % Struct containing the relevant information for the preconditioner.
    PS.droptol = droptol;
    PS.no_prec = false;
    PS.iter_refinement = false;
    PS.iter_refinement_maxit = 1;
    PS.stability_regularizer = 0;
    
    % Iterative refinement: external, internal or none.
    if ((NS.IR == "inner" || NS.IR == "inner_outer") && (min(NS.delta,NS.rho) < 5e-8)) 
        PS.iter_refinement = true;
        if (set_proximal_inner_IR)
            PS.iter_refinement_maxit = 1;
        end
        PS.stability_regularizer = 5e-8;
        stability_regularizer = PS.stability_regularizer;
    elseif (NS.prec_approach == "LDL_based_preconditioner" && ...
           (NS.solver == "minres" && NS.minres_setting == "augmented_system") && ...
                                             (min(NS.delta,NS.rho) < 5e-8))
        PS.stability_regularizer = 5e-8;
        stability_regularizer = PS.stability_regularizer;
    else
        stability_regularizer = NS.stability_regularizer;
    end      
    threshold = min(min(droptol,mu*droptol),1e0);
    %threshold = 0;
    n = NS.n;       m = NS.m;
    PS.n = n;       PS.m = m;
    PS.instability = false;
    PS.double_approximation = double_approximation;
    % =================================================================================================== %
    % building matrix E, which approximates Theta, and with it, build the approximate NE matrix M.
    % --------------------------------------------------------------------------------------------------- %
   
    B = true(n,1);
 
    E = 1./(NS.Theta_inv + NS.Q_struct.diag + stability_regularizer);
    PS.Q_barInv = E;
    N = (E<threshold);
%     out_inf = max(NS.d_inf,NS.p_inf);
%     if ((nnz(N) < 0.1*n) && (out_inf > 1e-4))
%         N = false(n,1);
%         if (out_inf > 5e-1)
%             force_drop = 0.25;
%         elseif (out_inf > 5e-3)
%             force_drop = 0.15;
%         else
%             force_drop = 0.1;
%         end
%         [~,I] = mink(E,ceil(force_drop*n));
%         N(I) = true;   
%     end
    B = xor(B,N);       % N = {1,...,n}\B.
    PS.nnz_B = nnz(B);
    if (NS.prec_approach == "Chol_based_preconditioner")
        B = and(B,C);
        PS.nnz_B = nnz(B);
        not_R = ~R;
        if ((~double_approximation) || (nnz(B) == 0) || (nnz(not_R) == 0)) % No additional approximation is needed. Proceed normally.
            E(N) = 0;
            PS.double_approximation = false;
            if (nnz(B) > 0)
                D = A(:,B);
                M = D*(spdiags(E(B),0,PS.nnz_B,PS.nnz_B)*D') + (NS.delta + stability_regularizer).*speye(m);
            else
                M = speye(m);
            end
            maxpiv = max(spdiags(M,0));
            if (~isnan(maxpiv) && ~isinf(maxpiv) )
                PS.instability = false;
                PS.maxpiv = maxpiv;
                if (PS.nnz_B > 0)
                    [PS.L_M,chol_flag,PS.P] = chol(M,'lower','vector');     % Cholesky factorization
                else
                    PS.L_M = speye(m);
                    chol_flag = 0;
                    PS.P = 1:m;
                    PS.no_prec = true;
                end
                PS.Pinv(PS.P) = 1:m;

                if (chol_flag ~= 0)
                    PS.instability = true;
                end
            else
                PS.instability = true;
                return;
            end
        else    % Double approximation is requested. Further approximate this preconditioner.
            [PS.R_NE_PS] = Reg_NE_preconditioner(A,E,A_tr,NS.delta + stability_regularizer,R,B);
            if (PS.R_NE_PS.instability_11 || PS.R_NE_PS.instability_22)
                PS.instability = true;
                return;
            end
        end
    elseif (NS.prec_approach == "LDL_based_preconditioner")
        PS.minres_setting = NS.minres_setting;
        if (NS.solver == "minres" && NS.minres_setting == "augmented_system")
            PS.iter_refinement_maxit = 1;       % No IR is allowed in this case.
            A_tilde = A;    nnz_N = nnz(N);
            pivot_threshold = 1e-10;
            if (nnz_N)
                A_tilde(:,N) = sparse(m,nnz_N);
                Q_tilde = sparse(n,n);
                Q_tilde(B,B) = Q(B,B);
                Q_tilde(N,N) = spdiags(NS.Q_struct.diag(N),0,nnz_N,nnz_N);
            else
                Q_tilde = Q;
            end
             K = [-Q_tilde - spdiags(NS.Theta_inv+stability_regularizer,0,n,n)               A_tilde';
                  A_tilde                 spdiags((NS.delta+stability_regularizer).*ones(m,1),0,m,m)];
        else
            pivot_threshold = 1e-10;
            A_tilde = A(:,B);
            Q_tilde = Q(B,B) + spdiags(NS.Theta_inv(B) + stability_regularizer,0,PS.nnz_B,PS.nnz_B);
            if (NS.solver == "minres")
                nnz_N = nnz(N);
                Q_N_tilde = NS.Q_struct.diag(N) + NS.Theta_inv(N) + stability_regularizer;
                Q_BN_tilde = sparse(n,n);
                Q_BN_tilde(N,N) = spdiags(Q_N_tilde,0,nnz_N,nnz_N);
                Q_BN_tilde(B,B) = Q_tilde;
                if (PS.nnz_B)
                    [PS.L_Q,chol_flag,PS.P_Q] = chol(Q_BN_tilde,'lower','vector');     % Cholesky factorization 
                    if (chol_flag ~= 0)
                        PS.instability = true;
                    end
                else
                    PS.L_Q = spdiags(Q_N_tilde.^(1/2),0,n,n);
                    PS.P_Q = 1:n;
                end
                PS.P_Qinv(PS.P_Q) = 1:n;
            end
            K = [ (-Q_tilde)                                                        A_tilde';
                  A_tilde           spdiags((NS.delta+stability_regularizer).*ones(m,1),0,m,m)];
          
        end
        pivots = (abs(spdiags(K,0)));
        maxpivot = max(pivots);
        minpivot = min(pivots);
        if (~isnan(minpivot) && (min(minpivot) > 0) && ~isinf(maxpivot))
            PS.maxpivot = maxpivot; PS.minpivot = minpivot;
            PS.instability = false;
            if (PS.nnz_B > 0)
                [PS.L_M,PS.D_M,PS.Pinv] = ldl(K,pivot_threshold,'vector');     % LDL' factorization
            else
                if (NS.solver == "minres" && NS.minres_setting == "augmented_system")
                    PS.L_M = speye(m+n);
                    PS.D_M = speye(m+n);
                    PS.Pinv = 1:(m+n);
                else
                    PS.L_M = speye(m);
                    PS.D_M = speye(m);
                    PS.Pinv = 1:m;
                end
                PS.no_prec = true;
            end
            if (NS.solver == "minres" && NS.minres_setting == "augmented_system")
                PS.P(PS.Pinv) = 1:(m+n);
            else
                PS.P(PS.Pinv) = 1:(m+PS.nnz_B);
            end
        else
            PS.instability = true;
            return;            
        end
    end
    % ____________________________________________________________________________________________________ %    
end
