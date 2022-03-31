function [NS] = Newton_matrix_data(iter,mu,A_struct,Q_struct,x,z_l,z_u,lb,ub,lb_vars,ub_vars,delta,rho,m,n,pivot_thr,solver,prec_approach,IR,p_inf,d_inf)
% ===================================================================================================================== %
% Newton_matrix_setting: Store all relevant information about the Newton matrix.
% --------------------------------------------------------------------------------------------------------------------- %
% NS = Newton_matrix_data(iter,x,z_l,z_u,lb,ub,lb_vars,ub_vars,delta,rho,m,n,pivot_thr,solver,approach,IR) returns a 
%      MATLAB struct that holds the relevant information of the Newton matrix, required for solving the step equations 
%      in the IPM.
% Author: Spyridon Pougkakiotis.
% _____________________________________________________________________________________________________________________ %
    NS = struct(); 
    NS.A_struct = A_struct;    NS.Q_struct = Q_struct;
    NS.m = m;                  NS.n = n;
    NS.x = x;                  NS.z_l = z_l;              NS.z_u = z_u;
    NS.lb = lb;                NS.ub = ub;
    NS.lb_vars = lb_vars;      NS.ub_vars = ub_vars;      NS.mu = mu;
    NS.Theta_inv = rho.*ones(n,1);     % lb_vars and ub_vars are mutually exclusive by construction.
                                       % It suffices to use a single Theta_inv matrix.
    NS.Theta_inv(lb_vars) = NS.Theta_inv(lb_vars) + z_l(lb_vars)./(x(lb_vars)-lb(lb_vars));
    NS.Theta_inv(ub_vars) = NS.Theta_inv(ub_vars) + z_u(ub_vars)./(ub(ub_vars)-x(ub_vars));
    NS.delta = delta;          NS.IPiter = iter;    NS.rho = rho;
    NS.pivot_thr = pivot_thr;   NS.prec_approach = prec_approach;
    NS.solver = solver; NS.p_inf = p_inf; NS.d_inf = d_inf;
    if (((IR == "outer") || (IR == "inner_outer")) && (min(delta,rho) < pivot_thr))
        NS.iter_refinement = true;
        NS.iter_refinement_maxit = 2;
        NS.stability_regularizer = 1e-10;  % Set if one needs additional stability (can hinder convergence).
        NS.IR_residual = false;
    else
        NS.iter_refinement = false;
        NS.IR_residual = false;
        NS.iter_refinement_maxit = 1;
        if (Q_struct.nnz)
            NS.stability_regularizer = 1e-10;
        else
            NS.stability_regularizer = 0;
        end
    end
    NS.IR = IR; 
end
% ******************************************************************************************************************** %
% END OF FILE.
% ******************************************************************************************************************** %
