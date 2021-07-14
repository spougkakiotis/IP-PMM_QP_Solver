function [w] = AS_multiplier(x,NS)
% ===================================================================================================================== %
% Augmented System Operator:
% --------------------------------------------------------------------------------------------------------------------- %
% [x] = AS_multiplier(x,NS,A_struct,Q_struct) takes as an input the struct containing the Newton blocks, as well as 
% a vector of size m, and returns the matrix-vector product of the augmented system's matrix by this vector.
% _____________________________________________________________________________________________________________________ %
    stability_regularizer = 0;
    if (NS.iter_refinement && ~(NS.IR_residual))
        stability_regularizer = NS.stability_regularizer;
    end
    w = zeros(NS.n+NS.m,1);
    w(1:NS.n) = -(NS.Theta_inv+stability_regularizer).*x(1:NS.n) - NS.Q_struct.Q(x(1:NS.n)) + NS.A_struct.A_tr(x(NS.n+1:NS.n+NS.m));
    w(NS.n+1:NS.n+NS.m) = NS.A_struct.A(x(1:NS.n)) + (NS.delta+stability_regularizer).*x(NS.n+1:NS.n+NS.m);
end

