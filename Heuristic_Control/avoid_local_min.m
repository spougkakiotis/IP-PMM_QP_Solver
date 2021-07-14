function [ctr_struct,reg_limit] = avoid_local_min(ctr_struct,rho,delta,reg_limit)
% ================================================================================================================ %
% avoid_local_min(ctr_struct,rho,delta,reg_limit): 
% ---------------------------------------------------------------------------------------------------------------- %
% This function checks whether the algorithm approaches dangerously close to a local minimum,
% causing the method to slow down. If so, it decreases the regularization threshold to 
% help the method escape it.
% ________________________________________________________________________________________________________________ %
        if (ctr_struct.no_primal_update > 5 && rho == reg_limit)
            reg_limit = 5*1e-13;
            ctr_struct.no_primal_update = 0;
            ctr_struct.no_dual_update = 0;
        elseif (ctr_struct.no_dual_update > 5 && delta == reg_limit)
            reg_limit = 5*1e-13;
            ctr_struct.no_primal_update = 0;
            ctr_struct.no_dual_update = 0;
        end
end