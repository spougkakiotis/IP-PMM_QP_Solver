function [alpha_x,alpha_z] = step_to_the_boundary(x,z_l,z_u,dx,dz_l,dz_u,lb,ub,lb_vars,ub_vars)
% ================================================================================================================ %
% step_to_the_boundary(x,z_l,z_u,dx,dz_l,dz_u,lb,ub,lb_vars,ub_vars): 
% ---------------------------------------------------------------------------------------------------------------- %
% This function computes the maximum step-length so that the new primal variable remains strictly within its 
% respective lower and upper bounds, while the dual slacks remain within the positive orthant.
% ________________________________________________________________________________________________________________ %
    % ============================================================================================================ %
    % Step in the lower-or-upper limit (primal) and non-negativity orthant (dual).
    % ------------------------------------------------------------------------------------------------------------ %
    tau = 0.995;        % tau in (0,1), to avoid taking the full Newton step and hitting the boundary.
    n = size(x,1);
    pos_idx = false(n,1);   neg_idx = false(n,1);
    neg_idz_l = false(n,1); neg_idz_u = false(n,1);
    pos_idx(ub_vars) = dx(ub_vars) > 0;
    neg_idx(lb_vars) = dx(lb_vars) < 0;
    neg_idz_l(lb_vars) = dz_l(lb_vars) < 0;
    neg_idz_u(ub_vars) = dz_u(ub_vars) < 0;
    alphamax_x = min([1;(lb(neg_idx)-x(neg_idx))./dx(neg_idx);...
                     (ub(pos_idx)-x(pos_idx))./dx(pos_idx)]);
    alphamax_z = min([1;-z_l(neg_idz_l)./dz_l(neg_idz_l);...
                     -z_u(neg_idz_u)./dz_u(neg_idz_u)]);
    alpha_x = tau*alphamax_x;
    alpha_z = tau*alphamax_z;
    % ____________________________________________________________________________________________________________ %
end