function [x,z_l,z_u,mu] = boundary_control(x,z_l,z_u,lb,ub, lb_vars,ub_vars,num_lb_vars,num_ub_vars,mu_prev)
% ============================================================================================= %
% boundary_control(x,z_l,z_u,lb_vars,ub_vars,num_lb_vars,num_ub_vars)
% --------------------------------------------------------------------------------------------- %
% This function takes as an input the newly computed IP-PMM iterate and checks whether 
% any of the constrained variables is dangerously close to the boundary. If so 
% is shifts it appropriately.
% _____________________________________________________________________________________________ %
    activate = false;
    if (num_lb_vars)
        if (min(x(lb_vars)-lb(lb_vars)) < 10*eps)
            x(lb_vars) = x(lb_vars) + 1e-11;
            activate = true;
        end
        if (min(z_l(lb_vars)) < 10*eps)
            z_l(lb_vars) = z_l(lb_vars) + 1e-11;
            activate = true;
        end
    end
    if (num_ub_vars)
        if (min(ub(ub_vars) - x(ub_vars)) < 10*eps)
            x(ub_vars) = x(ub_vars) - 1e-11;
            activate = true;
        end
        if (min(z_u(ub_vars)) < 10*eps)
            z_u(ub_vars) = z_u(ub_vars) + 1e-11;
            activate = true;
        end
    end
    if (activate)
        mu = (x(lb_vars)-lb(lb_vars))'*z_l(lb_vars) + (ub(ub_vars)-x(ub_vars))'*z_u(ub_vars);
       % mu = max(mu/(num_lb_vars+num_ub_vars),0.1*mu_prev);
        mu = mu/(num_lb_vars + num_ub_vars);
    else
        mu = mu_prev;
    end
end

