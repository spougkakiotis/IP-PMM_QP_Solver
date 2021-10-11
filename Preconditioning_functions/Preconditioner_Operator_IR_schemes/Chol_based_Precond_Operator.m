function [w] = Chol_based_Precond_Operator(x,PS,solver)
% ===================================================================================================================== %
% Preconditioner Operator:
% --------------------------------------------------------------------------------------------------------------------- %
% [w] = Chol_based_Precond_Operator(x,PS,solver) takes as an input the struct containing the preconditioner blocks, as well as 
% a vector of size n or n+m, and returns the matrix-vector product of the inverse preconditioner by this vector. 
% The latter can be approximated via an Iterative Refinement approach.
% Two Cholesky-based preconditioners are supported:
%       -> Normal equations approximation + pcg
%       -> Augmented system approximation + minres (via block diagonal preconditioning), with \tilde{Q} = diag(Q).
% _____________________________________________________________________________________________________________________ % 
    for i = 1:PS.iter_refinement_maxit
        if (i == 1)
            rhs = x;
        end
        if (solver == "pcg")
            if (i > 1)
                rhs = x + PS.stability_regularizer.*w;             
            end
            if (~PS.double_approximation || PS.no_prec)
                w = Preconditioner_backsolve(PS.L_M,PS.P,PS.Pinv,rhs);
            else
                w = RNE_precond_operator(rhs,PS.R_NE_PS);
            end
        elseif (solver == "minres")
            if (i == 1)
                w = zeros(PS.n+PS.m,1);
            else
                rhs = x + PS.stability_regularizer.*w;
            end
            w(1:PS.n) = PS.Q_barInv .* rhs(1:PS.n);
            if (~PS.double_approximation || PS.no_prec)
                w(PS.n+1:PS.n+PS.m,1) = Preconditioner_backsolve(PS.L_M,PS.P,PS.Pinv,rhs(PS.n+1:PS.n+PS.m,1));
            else
                w(PS.n+1:PS.n+PS.m,1) = RNE_precond_operator(rhs(PS.n+1:PS.n+PS.m,1),PS.R_NE_PS);
            end
        else
            error('Incorrect Input argument: solver.');
        end
    end
end

