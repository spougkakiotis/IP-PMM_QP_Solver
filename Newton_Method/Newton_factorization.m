 function NS = Newton_factorization(A,A_tr,Q,x,z_l,z_u,delta,rho,lb,ub,lb_vars,ub_vars,pivot_thr,IR)
% ==================================================================================================================== %
% Newton_factorization: Factorize the Newton matrix
% -------------------------------------------------------------------------------------------------------------------- %
% NS = Newton_factorization(A,A_tr,Q,x,z,delta,rho,pos_vars,free_vars) returns a MATLAB struct that holds the
%      factorization of the Newton matrix for solving the step equations in
%      the IPM, as well as relevant information concerning failure.
% Factorization Method
  % --------------------
  % 1: augmented system, LDL' factorization.
% 
% Author: Spyridon Pougkakiotis.
% ==================================================================================================================== %
[m, n] = size(A);
NS = struct();
% ==================================================================================================================== %
% LDL' factorization of KKT matrix
% -------------------------------------------------------------------------------------------------------------------- %
% Perform the same reduction as above but calculate its symmetric
% indefinite factorization
%
%   K(pp,pp) = L*D*L'.
%
% MATLAB uses MA57 when K is sparse, which is not available in OCTAVE. 
% -------------------------------------------------------------------------------------------------------------------- %
NS.x = x;                  NS.z_l = z_l;              NS.z_u = z_u;
NS.lb = lb;                NS.ub = ub;
NS.lb_vars = lb_vars;      NS.ub_vars = ub_vars;
if ((IR == "outer") && (min(delta,rho) < pivot_thr))
    NS.iter_refinement = true;
    NS.iter_refinement_maxit = 2;
    NS.stability_regularizer = pivot_thr*(5e0);
    NS.IR_residual = false;
else
    NS.iter_refinement = false;
    NS.IR_residual = false;
    NS.iter_refinement_maxit = 1;
    NS.stability_regularizer = 0;
end
NS.IR = IR; 
Theta_inv_l = zeros(n,1);  Theta_inv_u = zeros(n,1);
Theta_inv_l(lb_vars) = z_l(lb_vars)./(x(lb_vars)-lb(lb_vars));
Theta_inv_u(ub_vars) = z_u(ub_vars)./(ub(ub_vars)-x(ub_vars));
K = [-Q-spdiags(Theta_inv_l + Theta_inv_u + (rho+NS.stability_regularizer).*ones(n,1),0,n,n),                   A_tr; 
      A,                                                  spdiags((delta+NS.stability_regularizer).*ones(m,1),0,m,m)]; 
[NS.L,NS.D,NS.pp] = ldl(K,pivot_thr,'vector'); %Small pivots allowed, to avoid 2x2 pivots.
% ==================================================================================================================== %  
 
% ******************************************************************************************************************** %
% END OF FILE.
% ******************************************************************************************************************** %
end
