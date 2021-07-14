function [dx,dy,dz_l,dz_u,instability,Krylov_iter,drop_direction,Reinforce_Inner_IR] = Newton_itersolve(fid,pred,NS,PS,res_p,res_d,res_mu_l,res_mu_u,maxit,tol,printlevel)
% ======================================================================================================================== %
% Newton_itersolve    Solve linear system with a Krylov method, by using a given preconditioner.
% ------------------------------------------------------------------------------------------------------------------------ %
% OUTPUT:
%  [dx,dy,dz,instability,iter,drop_direction,Reinforce_Inner_IR] = Newton_itersolve(fid,pred,NS,PS,res_p,res_d,res_mu,maxit,tol,printlevel)
%  Returns the Newton direction, a parameter indicating critical ill-conditioning, the number of Krylov iterations,
%  as well as a request to drop the direction if accuracy is far from being reached.
%
% Author: Spyridon Pougkakiotis.
% ________________________________________________________________________________________________________________________ %
    m = size(res_p,1);
    n = size(res_d,1);
    instability = false;        drop_direction = false;         Reinforce_Inner_IR = false;
    dx = zeros(n,1);            dy = zeros(m,1);
    dz_l = zeros(n,1);          dz_u = zeros(n,1);
    temp_res_l = zeros(n,1);    temp_res_u = zeros(n,1);
    out_inf = min(NS.d_inf,NS.p_inf);
    % ==================================================================================================================== %
    % Solve KKT system with LDL' factors.
    % -------------------------------------------------------------------------------------------------------------------- %
    temp_res_l(NS.lb_vars) =  res_mu_l(NS.lb_vars)./(NS.x(NS.lb_vars)-NS.lb(NS.lb_vars));
    temp_res_u(NS.ub_vars) =  res_mu_u(NS.ub_vars)./(NS.ub(NS.ub_vars)-NS.x(NS.ub_vars));
    rhs = [res_d-temp_res_l+temp_res_u; res_p];
    solver = NS.solver;
    if (solver == "pcg")
        Theta = 1./(NS.Theta_inv+NS.Q_struct.diag + NS.stability_regularizer);
        rhs_y = NS.A_struct.A(Theta.*rhs(1:n,1)) + rhs(n+1:n+m,1);
    else
        Theta = 1./(NS.Theta_inv);
    end
    % ____________________________________________________________________________________________________________________ %
    
    % ==================================================================================================================== %
	% Call the respective solver the problem using the constructed  preconditioner, whose parts are stored in struct PS.
	% -------------------------------------------------------------------------------------------------------------------- %
    warn_stat = warning;
    warning('off','all');
    Krylov_iter = 0;

    if (solver == "pcg")
        accuracy_bound = 1e-4;
        if (NS.prec_approach == "LDL_based_preconditioner")
            solver = NS.prec_approach;
        end
        for i = 1:NS.iter_refinement_maxit
            if (i == 1)
                res_rhs = rhs_y;
                IR_residual = rhs;
                lhs = zeros(n+m,1);
            else
                NS.IR_residual = true;                                                      % Turns off stability.
                IR_residual = rhs - AS_multiplier(lhs,NS);                                  % Computes the residual.
                NS.IR_residual = false;
                if (norm(IR_residual) <= 1e-8)                                              % No need to refine.
                    break;
                end
                res_rhs = IR_residual(n+1:n+m) + NS.A_struct.A(Theta.*IR_residual(1:n));    % PCG right-hand-side.
            end
            tol = max(1e-10,tol/max(1,norm(res_rhs,'Inf')));
            accuracy_bound = max(tol*1e1, accuracy_bound);
            [lhs_y, flag, res, iter] = pcg(@(x) NE_multiplier(x,NS), res_rhs, tol, maxit, @(x) Precond_Operator(x,PS,solver));
            lhs_x = (Theta).*(-IR_residual(1:n) + NS.A_struct.A_tr(lhs_y));
            lhs = lhs + [lhs_x; lhs_y];
            Krylov_iter = Krylov_iter + iter;
            if ((res > accuracy_bound) || (flag >= 2 && flag ~= 3))
                break;
            elseif (flag == 3 && PS.iter_refinement)
                Reinforce_Inner_IR = true;   % If we use stability on preconditioner, it must be reinforced.
            end
        end
    elseif (solver == "minres")
        accuracy_bound = 1e-3;
        if (NS.prec_approach == "LDL_based_preconditioner")
            solver = NS.prec_approach;
        end
        for i = 1:NS.iter_refinement_maxit
            if (i == 1)
                IR_residual = rhs;
                lhs = zeros(n+m,1);
            else
                NS.IR_residual = true;
                IR_residual = rhs - AS_multiplier(lhs,NS);
                NS.IR_residual = false;
                if (norm(IR_residual) <= 1e-8)
                    break;
                end
            end
            tol = max(1e-10,(tol)/max(1,norm(IR_residual,'Inf')));
            %tol = 0.5*min(1e-3,out_inf);
            %tol = max(tol,1e-4);
            accuracy_bound = max(tol*1e1, accuracy_bound);
            [lhs_c, flag, res, iter] = minres(@(x) AS_multiplier(x,NS), IR_residual, tol, maxit, @(x) Precond_Operator(x,PS,solver));
            lhs = lhs + lhs_c;
            Krylov_iter = Krylov_iter + iter;
            if ((res > accuracy_bound) || (flag >= 2 && flag ~= 3))
                break;
            elseif (flag == 3 && PS.iter_refinement)
                Reinforce_Inner_IR = true;   % If we use stability on preconditioner, it must be reinforced.
            end
        end
    end
    warning(warn_stat);
    % ____________________________________________________________________________________________________________________ %
    % ==================================================================================================================== %
    % Compute the Newton directions and report the relevant statistics.
    % -------------------------------------------------------------------------------------------------------------------- %
    if (flag > 0) % Something went wrong, so we assume that the preconditioner is not good enough -> increase quality.
        Krylov_iter = max(maxit,Krylov_iter);
        if (flag >= 2 && flag ~= 3 && flag ~= 5)
            instability = true;
            fprintf('Instability detected during the iterative method. flag = %d.\n',flag);
            return;
        elseif ((res > accuracy_bound))
            drop_direction = true;
            return;
        end
    end
    
    if (printlevel >= 3)
        if (~pred)
            fprintf(fid,'                               ***Krylov method:Corrector***                                           \n');
            fprintf(fid,'Krylov Iterations                Krylov Flag                 Residual               Droptol            \n');
            fprintf(fid,'%4d                           %4d                         %9.2e              %9.2e                 \n',Krylov_iter,flag,res,PS.droptol);
            fprintf(fid,'_______________________________________________________________________________________________________\n');
        else
            fprintf(fid,'-------------------------------***Krylov method:Predictor***-------------------------------------------\n');
            fprintf(fid,'Krylov Iterations                Krylov Flag                 Residual               Droptol            \n');
            fprintf(fid,'%4d                           %4d                         %9.2e              %9.2e                 \n',Krylov_iter,flag,res,PS.droptol);
        end
    end
   
    if (nnz(isnan(lhs)) > 0 || nnz(isinf(lhs)) > 0 || (max(lhs) == 0 && min(lhs) == 0)) % Check for ill-conditioning.
        instability = true;
        Krylov_iter = max(maxit,Krylov_iter);
        fprintf('Instability detected during the iterative method.\n');
        return;
    end
    dx = lhs(1:n,1);
    dy = lhs(n+1:n+m,1);
    dz_l(NS.lb_vars) = (res_mu_l(NS.lb_vars) - NS.z_l(NS.lb_vars).*dx(NS.lb_vars))./(NS.x(NS.lb_vars)-NS.lb(NS.lb_vars));
    dz_u(NS.ub_vars) = (res_mu_u(NS.ub_vars) + NS.z_u(NS.ub_vars).*dx(NS.ub_vars))./(NS.ub(NS.ub_vars)-NS.x(NS.ub_vars));
    % ____________________________________________________________________________________________________________________ %
end 