function [dx,dy,dz_l,dz_u,instability] = Newton_backsolve(NS,res_p,res_d,res_mu_l,res_mu_u)
% ==================================================================================================================== %
% Newton_backsolve    Solve linear system with factorized matrix, by using backward substitution.
% -------------------------------------------------------------------------------------------------------------------- %
% OUTPUT:
%  [dx,dy,dz,instability] = newtonsolve(NS,res_p,res_d,res_mu,A,A_tr,pos_vars,free_vars)
%  i.e. the Newton direction and a boolean parameter indicating critical ill-conditioning.
%
% Author: Spyridon Pougkakiotis.
% ____________________________________________________________________________________________________________________ %
m = size(res_p,1);
n = size(res_d,1);
instability = false;
dx = zeros(n,1);            dy = zeros(m,1);
dz_l = zeros(n,1);          dz_u = zeros(n,1);
temp_res_l = zeros(n,1);    temp_res_u = zeros(n,1);
% ==================================================================================================================== %
% Solve KKT system with LDL' factors.
% -------------------------------------------------------------------------------------------------------------------- %
temp_res_l(NS.lb_vars) =  res_mu_l(NS.lb_vars)./(NS.x(NS.lb_vars)-NS.lb(NS.lb_vars));
temp_res_u(NS.ub_vars) =  res_mu_u(NS.ub_vars)./(NS.ub(NS.ub_vars)-NS.x(NS.ub_vars));
rhs = [res_d-temp_res_l+temp_res_u; res_p];
for i = 1:NS.iter_refinement_maxit
    if (i == 1)
        IR_residual = rhs;
        lhs_prev = zeros(n+m,1);
    else
        lhs_prev = lhs;
        non_regularized_product = lhs_prev(NS.pp);
        non_regularized_product = (NS.L*(NS.D*(NS.L'*non_regularized_product))) ...
           - [-NS.stability_regularizer.*ones(n,1); NS.stability_regularizer.*ones(m,1)].*non_regularized_product;
        non_regularized_product(NS.pp) = non_regularized_product;
        IR_residual = rhs - non_regularized_product;                                % Computes the residual.
        if (norm(IR_residual) <= 1e-8)                                              % No need to refine.
            break;
        end
    end
    warn_stat = warning;
    warning('off','all');
    lhs = NS.L'\(NS.D\(NS.L\IR_residual(NS.pp)));
    if (nnz(isnan(lhs)) > 0 || nnz(isinf(lhs)) > 0)
        instability = true;
        return;
    end
    warning(warn_stat);
    lhs(NS.pp) = lhs;
    lhs = lhs_prev + lhs;        
end
dx = lhs(1:n,1);
dy = lhs(n+1:n+m,1);
dz_l(NS.lb_vars) = (res_mu_l(NS.lb_vars) - NS.z_l(NS.lb_vars).*dx(NS.lb_vars))./(NS.x(NS.lb_vars)-NS.lb(NS.lb_vars));
dz_u(NS.ub_vars) = (res_mu_u(NS.ub_vars) + NS.z_u(NS.ub_vars).*dx(NS.ub_vars))./(NS.ub(NS.ub_vars)-NS.x(NS.ub_vars));
% ____________________________________________________________________________________________________________________ %

% ******************************************************************************************************************** %
% END OF FILE.
% ******************************************************************************************************************** %
end 