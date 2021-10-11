function [w] = Precond_Operator(x,PS,solver)
% ===================================================================================================================== %
% Preconditioner Operator:
% --------------------------------------------------------------------------------------------------------------------- %
% [w] = Precond_Operator(x,PS,solver) takes as an input the struct containing the preconditioner blocks, as well as 
% a vector of size n+m, and returns the matrix-vector product of the inverse preconditioner by this vector.
% _____________________________________________________________________________________________________________________ % 
    if (solver == "LDL_based_preconditioner")
        [w] =  LDL_based_Precond_Operator(x,PS);
    elseif (solver == "pcg" || solver == "minres")
        [w] = Chol_based_Precond_Operator(x,PS,solver);
    else
        error('Incorrect Input argument: solver.')    
    end
end

