function [x] = NE_multiplier(x,NS)
% ==================================================================================================================== %
% Normal Equations Operator:
% --------------------------------------------------------------------------------------------------------------------- %
% [x] = NE_multiplier(x,NS,A_struct,Q_struct) takes as an input the struct containing the Newton blocks, as well as 
% a vector of size m, and returns the matrix-vector product of the normal equations' matrix by this vector.
% _____________________________________________________________________________________________________________________ %
    stability_regularizer = 0;
    if (NS.iter_refinement && ~(NS.IR_residual))
        stability_regularizer = NS.stability_regularizer;
    end
    w = NS.A_struct.A_tr(x);
    w = (1./(NS.Theta_inv+NS.Q_struct.diag + stability_regularizer)).*w;  
    w = NS.A_struct.A(w);
    x = w + (NS.delta+stability_regularizer).*x;
end

