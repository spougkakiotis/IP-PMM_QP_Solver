function [x,y,z_l,z_u] = Mehrotra_Warm_Start(A_struct, Q_struct, b, c, lb, ub, lb_vars, ub_vars)
% ==================================================================================================================== %
% Initialization - Mehrotra-like Initial Point for QP:
% -------------------------------------------------------------------------------------------------------------------- %
% Choose an initial starting point (x,y,z_l,z_u). For that, we ignore the box constraints, 
% and solve the relaxed regularized optimization problem (which has a closed form solution). Then,
% we shift the solution, to respect the box constraints. The point is expected to be well centered.
%
% Author: Spyridon Pougkakiotis, January 2021, Edinburgh.
% ____________________________________________________________________________________________________________________ %
    % ================================================================================================================ %
    % Use PCG to solve two least-squares problems for efficiency (along with the Jacobi preconditioner). 
    % ---------------------------------------------------------------------------------------------------------------- %
    m = size(b,1); n = size(c,1);
    delta_0 = 1e-3;
    D = A_struct.NE_diag + delta_0;
    Jacobi_Prec = @(x) (1./D).*x;
    NE_fun = @(x) (A_struct.A(A_struct.A_tr(x)) + delta_0.*x);
    if (norm(b) <= 1e-6) b = b + 1e-3.*ones(m,1); end
    x = pcg(NE_fun,b,1e-3,min(200,m),Jacobi_Prec);
    x = A_struct.A_tr(x);
    y = pcg(NE_fun,A_struct.A(c+Q_struct.Q(x)),1e-3,min(200,m),Jacobi_Prec);
    z = (c + Q_struct.Q(x) - A_struct.A_tr(y));
    z_l = zeros(n,1); z_u = zeros(n,1);
    z_l(lb_vars) = z(lb_vars); z_u(ub_vars) = -z(ub_vars);   %lb_vars and ub_vars are mutually exclusive (IP-PMM format).
    % ________________________________________________________________________________________________________________ %
    num_lb_vars = nnz(lb_vars); num_ub_vars = nnz(ub_vars);
    % ================================================================================================================ %
    % Avoid solutions too close to the boundary.
    % ---------------------------------------------------------------------------------------------------------------- %
    if (norm(x(lb_vars)-lb(lb_vars)) <= 10^(-4)) 
        x(lb_vars) = 0.1.*ones(num_lb_vars,1) + lb(lb_vars); % 0.1 is chosen arbitrarily
    end
    if (norm(ub(ub_vars)-x(ub_vars)) <= 10^(-4)) 
        x(ub_vars) = -0.1.*ones(num_ub_vars,1) + ub(ub_vars); % 0.1 is chosen arbitrarily
    end
    if (norm(z_l(lb_vars)) <= 10^(-4))
        z_l(lb_vars) = 0.1.*ones(num_lb_vars,1); % 0.1 is chosen arbitrarily
    end
    if (norm(z_u(ub_vars)) <= 10^(-4))
        z_u(ub_vars) = 0.1.*ones(num_ub_vars,1); % 0.1 is chosen arbitrarily
    end
    % ________________________________________________________________________________________________________________ %
    
    
    % ================================================================================================================ %
    % Shift the solution to obtain something within the box bounds. Distinguish three cases for x:
    %           1. Has only a lower bound -> Shift-up according to Mehrotra's heuristic,
    %           2. Has only an upper bound -> Shift-down according to Mehrotra's heuristic,
    %           3. Has upper and lower bounds -> Take the average value of the bounds.
    % ---------------------------------------------------------------------------------------------------------------- %
    delta_x_l = max(-1.5*min(x(lb_vars)-lb(lb_vars)),0);
    delta_x_u = max(-1.5*min(ub(ub_vars)-x(ub_vars)),0);
    delta_z_l = max(-1.5*min(z_l(lb_vars)), 0);
    delta_z_u = max(-1.5*min(z_u(ub_vars)), 0);
    temp_product_l = (x(lb_vars)-lb(lb_vars) + delta_x_l.*ones(num_lb_vars,1))'*(z_l(lb_vars) ...
                   + delta_z_l.*ones(num_lb_vars,1));
    temp_product_u = (ub(ub_vars)-x(ub_vars) + ... 
                      delta_x_u.*ones(num_ub_vars,1))'*(z_u(ub_vars)  + delta_z_u.*ones(num_ub_vars,1));

    delta_x_l_bar = delta_x_l + (0.5*temp_product_l)/(sum(z_l(lb_vars),1) + num_lb_vars*delta_z_l);
    delta_x_u_bar = delta_x_u + (0.5*temp_product_u)/(sum(z_u(ub_vars),1) + num_ub_vars*delta_z_u);
    delta_z_l_bar = delta_z_l + (0.5*temp_product_l)/(sum(x(lb_vars)-lb(lb_vars),1)+num_lb_vars*delta_x_l);
    delta_z_u_bar = delta_z_u + (0.5*temp_product_u)/(sum(ub(ub_vars)-x(ub_vars),1)+num_ub_vars*delta_x_u);
        
    x(lb_vars) = x(lb_vars) + delta_x_l_bar.*ones(num_lb_vars,1);   % We have either lower OR upper bound.
    x(ub_vars) = x(ub_vars) - delta_x_u_bar.*ones(num_ub_vars,1);
    z_l(lb_vars) = z_l(lb_vars) + delta_z_l_bar.*ones(num_lb_vars,1);
    z_u(ub_vars) = z_u(ub_vars) + delta_z_u_bar.*ones(num_ub_vars,1);
    % ________________________________________________________________________________________________________________ %
    if (issparse(x))    x = full(x);     end
    if (issparse(z_l))  z_l = full(z_l); end
    if (issparse(z_u))  z_u = full(z_u); end
    if (issparse(y))    y = full(y);     end
end

