 function NS = Newton_factorization(A,A_tr,Q,x,z_l,z_u,delta,rho,lb,ub,lb_vars,ub_vars,pivot_thr)
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
Theta_inv_l = zeros(n,1);  Theta_inv_u = zeros(n,1);
Theta_inv_l(lb_vars) = z_l(lb_vars)./(x(lb_vars)-lb(lb_vars));
Theta_inv_u(ub_vars) = z_u(ub_vars)./(ub(ub_vars)-x(ub_vars));
K = [-Q-spdiags(Theta_inv_l + Theta_inv_u + rho.*ones(n,1),0,n,n),                           A_tr; 
      A,                                                           spdiags(delta.*ones(m,1),0,m,m)]; 
[NS.L,NS.D,NS.pp] = ldl(K,pivot_thr,'vector'); %Small pivots allowed, to avoid 2x2 pivots.
% ==================================================================================================================== %  
 
% ******************************************************************************************************************** %
% END OF FILE.
% ******************************************************************************************************************** %
end
