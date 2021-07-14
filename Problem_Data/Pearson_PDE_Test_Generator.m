function [solution_statistics] = Pearson_PDE_Test_Generator(problem_choice,tol,IPM_maxit,printlevel,la_mode,fid)
% ====================================================================================================== %
% This function takes an integer as an input, which specifies the problem 
% that is being solved by SSN-PMM. Then, it generates and return the relevant data 
% required to solve the problem, in the form of a structure.
% ------------------------------------------------------------------------------------------------------ %
% Input
%        problem_choice: 1 -> Poisson Control with L^2 regularization and no bounds,
%                        2 -> Poisson Control with L^2 regularization and control bounds,
%                        3 -> Convection Diffusion with H^1-regularization and state and control bounds,
%                        4 -> Poisson Control with L^1/L^2-regularization and control bounds,
%                        5 -> Convection Diffusion with L^1/L^2-regularization and state and control bounds,
% Output
%        solution_statistics: A structure of structures containing all relevant problem data.
% ______________________________________________________________________________________________________ %
    solution_statistics = struct();
    solution_statistics.total_iters = 0; solution_statistics.total_time = 0;
    solution_statistics.problem_converged = 0; solution_statistics.tol = tol;
    solution_statistics.total_Krylov_iters = 0;
    solution_statistics.objective_value = 0;    % To keep objective values.
    solution_statistics.status = 0;          % To keep convergence status.
    solution_statistics.max_nnzL = 0;            % To keep the maximum nnz of a Cholesky factor found.

    grid_type = 1; % uniform grid
    % ================================================================================================== %
    % Request input parameters from the user or use the default ones.
    % -------------------------------------------------------------------------------------------------- %
    fprintf('Should the default parameters be used?\n');
    default = input('Type 1 for default parameters, or anything else to manually include them.\n');
    if (default == 1)
            nc = 6;
            gamma = 1e-2;
    else 
        fprintf('Choose the number of uniform grid points.\n');
        while(true)
            nc = input('Type an integer k in [3,15] to enforce 2^k+1 grid points in each direction.\n');
            if (isinf(nc) || isnan(nc) || floor(nc)~= nc || nc >= 16 || nc <= 2)
                fprintf('Incorrect input argument. Please type an integer k between 3 and 15.\n');
            else
                break;
            end
        end
        fprintf('Choose a value for the smooth regularization parameter.\n');
        while(true)
            gamma = input('Type a double value in the form 1e-k, where k must be in [0,12].\n');
            if (isinf(gamma) || isnan(gamma) || ~isa(gamma,'double') ||...
                gamma > 1e-0 || gamma < 1e-12)
                fprintf('Incorrect input argument.\n');
            else
                break;
            end
        end
        if (problem_choice == 4 || problem_choice == 5)
            fprintf('Choose a value for the L1 regularization parameter\n');
            while(true)
                beta = input('Type a double value in the form 1e-k, where k must be in [0,12].\n');
                if (isinf(beta) || isnan(beta) || ~isa(beta,'double') ||...
                    beta > 1e-0 || beta < 1e-12)
                    fprintf('Incorrect input argument.\n');
                else
                    break;
                end
            end
        end
        if (problem_choice == 2 || problem_choice == 4)
            fprintf('Choose two uniform lower and upper bounds for the control.\n');
            while (true)
                u_alpha = input('Type the lower bound as a real number.\n');
                u_beta = input('Type the upper bound as a real number.\n');
                if ( (u_alpha < -1) || (u_alpha > 0) || (u_alpha > u_beta) || (u_beta < 0) || (u_beta > 1))
                    fprintf('Incorrect input argument.\n');
                else
                    break;
                end
            end
        end
        if (problem_choice == 3 || problem_choice == 5)
           fprintf('Choose two uniform lower and upper bounds for the state and the control.\n');
            while (true)
                u_alpha = input('Type the lower bound as a real number in [-2, 0].\n');
                u_beta = input('Type the upper bound as a real number in [0, 2].\n');
                if ( (u_alpha < -2) || (u_alpha > 0) || (u_alpha >= u_beta) || (u_beta < 0) || (u_beta > 2))
                    fprintf('Incorrect input argument.\n');
                    continue;
                end
                y_alpha = input('Type the lower bound as a real number in [-1, 0].\n');
                y_beta = input('Type the upper bound as a real number in [0, 1].\n');
                if ( (y_alpha < -1) || (y_alpha > 0) || (y_alpha >= y_beta) || (y_beta < 0) || (y_beta > 1))
                    fprintf('Incorrect input argument.\n');
                else
                    break;
                end
            end   
            fprintf('Choose a real value for the diffusion coefficient.\n');
            while (true)
                epsilon = input('Type the value as a real number in [0.01,0.05].\n');
                if (epsilon < 0.1 || epsilon > 0.5)
                    fprintf('Incorrect input argument.\n');
                else
                    break;
                end
            end
        end
    end
    % __________________________________________________________________________________________________ %
     fileID = fopen('./Output_files/Pearson_PDE_tabular_format_IP_PMM_runs.txt','a+');

     % How large is resulting stiffness or mass matrix?
     np = (2^nc+1)^2; % entire system is thus 3*np-by-3*np

     % Compute matrices specifying location of nodes
    [x_1,x_2,x_1x_2,bound,mv,mbound] = square_domain_x(nc,grid_type);

    if (problem_choice == 1)
        % ============================================================================================== %
        % Poisson Control with L^2 regularization and no bounds.
        % ---------------------------------------------------------------------------------------------- %
        problem_name = "Unbounded_Poisson_Control";
        % Compute connectivity, stiffness and mass matrices (D and J)
        [ev,ebound] = q1grid(x_1x_2,mv,bound,mbound);
        [D,J_y] = femq1_diff(x_1x_2,ev);
        O = sparse(np,np);
        
        % Specify vectors relating to desired state and Dirichlet BCs
        yhat_vec = 1-(1+4*gamma*pi^4)/(2*pi^2)*sin(pi*x_1x_2(:,1)).*sin(pi*x_1x_2(:,2));
        bc_nodes = ones(length(bound),1);

        % Initialize RHS vector corresponding to desired state
        Jyhat = J_y*yhat_vec;                   
        Jf = zeros(np,1);

        % Enforce Dirichlet BCs on state, and zero Dirichlet BCs on adjoint
        [D,b] = nonzerobc_input(D,Jf,x_1x_2,bound,bc_nodes);
        [J_y,Jyhat] = nonzerobc_input(J_y,Jyhat,x_1x_2,bound,bc_nodes);
        J_constr = J_y; for i = 1:length(bound), J_constr(bound(i),bound(i)) = 0; end

        
        % Use the produced data to construct the relevant structures that SSN-PMM will use.
        obj_const_term = (1/2).*(yhat_vec'*Jyhat);      % Constant term included in the objective.
        c = [-Jyhat; zeros(np,1)];
  
        A = [D J_constr];
        Q = [J_y O; O gamma.*J_y];
        lb = -Inf.*ones(2*np,1);
        ub = Inf.*ones(2*np,1);
  
        [solution_struct] = Set_Up_IP_PMM(A,Q,b,c,obj_const_term,lb,ub,tol,IPM_maxit,printlevel,la_mode,fid);
        % Compare computed solutions to exact solutions
        y_ex = 1-1/(2*pi^2)*sin(pi*x_1x_2(:,1)).*sin(pi*x_1x_2(:,2));
        u_ex = sin(pi*x_1x_2(:,1)).*sin(pi*x_1x_2(:,2));
        p_ex = -gamma*sin(pi*x_1x_2(:,1)).*sin(pi*x_1x_2(:,2));
        y_relerr = norm(solution_struct.x(1:np)-y_ex)/norm(y_ex);
        u_relerr = norm(solution_struct.x(np+1:2*np)-u_ex)/norm(u_ex);
        p_relerr = norm(-solution_struct.y-p_ex)/norm(p_ex);
        ex_obj_val = c'*[y_ex;u_ex] + (1/2)*([y_ex;u_ex]'*(Q*[y_ex;u_ex])) + obj_const_term;
        obj_gap = solution_struct.obj_val - ex_obj_val;
        fprintf(['The objective value is %.2e. The gap w.r.t. the optimal objective is:\n',...
                 'Objective Gap = %.2e.\n'],solution_struct.obj_val, obj_gap);

        fprintf(['The relative errors w.r.t. the optimal solution are: \n', ...
                '|y-y_ex|/|y_ex| = %.2e, |u-u_ex|/|u_ex| = %.2e, |v-v_ex|/|v_ex| = %.2e.\n'],...
                y_relerr, u_relerr, p_relerr);      
        solution_statistics.y_ex = y_ex;
        solution_statistics.u_ex = u_ex;
        solution_statistics.p_ex = p_ex;
        solution_statistics.ex_obj_val = ex_obj_val;

        % ______________________________________________________________________________________________ %
    elseif (problem_choice == 2)
        % ============================================================================================== %
        % Poisson Control with L^2 regularization and bounds on the control.
        %     For this problem, continuous optimality conditions are (L the Laplacian):
        %         -Ly+u=f, -Lp=yhat-y, u=proj_[-alpha,alpha](-p/gamma)
        %     These give the exact solutions as below.
        % ---------------------------------------------------------------------------------------------- %   
        problem_name = "Bounded_Poisson_Control";
        if (default == 1)
            u_alpha = -0.5; % control is constrained to be within range [-alpha,alpha]
            u_beta = 0.5;
        end

        % Compute connectivity, stiffness and mass matrices (D and J)
        [ev,ebound] = q1grid(x_1x_2,mv,bound,mbound);
        [D,J_y] = femq1_diff(x_1x_2,ev);
        O = sparse(np,np);
        
        
        
        % Specify vectors relating to desired state, source term and Dirichlet BCs
        yhat_vec = 1-(1+4*gamma*pi^4)/(2*pi^2)*sin(pi*x_1x_2(:,1)).*sin(pi*x_1x_2(:,2));
        f_vec = -sin(pi*x_1x_2(:,1)).*sin(pi*x_1x_2(:,2));
        f_vec = f_vec+max(min(-f_vec,u_beta),u_alpha);
        bc_nodes = ones(length(bound),1);

        % Initialize RHS vector corresponding to desired state
        Jyhat = J_y*yhat_vec;
        Jf = J_y*f_vec;

        % Enforce Dirichlet BCs on state, and zero Dirichlet BCs on adjoint
        [D,b] = nonzerobc_input(D,Jf,x_1x_2,bound,bc_nodes);
        [J_y,Jyhat] = nonzerobc_input(J_y,Jyhat,x_1x_2,bound,bc_nodes);
        J_constr = J_y; for i = 1:length(bound), J_constr(bound(i),bound(i)) = 0; end

        obj_const_term = (1/2).*(yhat_vec'*Jyhat);      % Constant term included in the objective.
        c = [-Jyhat; zeros(np,1)];
        
        A = [D J_constr];
        Q = [J_y O; O gamma.*J_y];
        lb = -Inf.*ones(2*np,1);
        ub = Inf.*ones(2*np,1);
        lb(np+1:2*np) = u_alpha.*ones(np,1);
        ub(np+1:2*np) = u_beta.*ones(np,1);

        
        [solution_struct] = Set_Up_IP_PMM(A,Q,b,c,obj_const_term,lb,ub,tol,IPM_maxit,printlevel,la_mode,fid);
 
        % Exact solutions for this problem
        y_ex = 1-1/(2*pi^2)*sin(pi*x_1x_2(:,1)).*sin(pi*x_1x_2(:,2));
        u_ex = max(min(sin(pi*x_1x_2(:,1)).*sin(pi*x_1x_2(:,2)),u_beta),u_alpha);
        p_ex = -gamma*sin(pi*x_1x_2(:,1)).*sin(pi*x_1x_2(:,2));

        
        y_relerr = norm(solution_struct.x(1:np)-y_ex)/norm(y_ex);
        u_relerr = norm(solution_struct.x(np+1:2*np)-u_ex)/norm(u_ex);
        p_relerr = norm(-solution_struct.y(1:np)-p_ex)/norm(p_ex);
        ex_obj_val = c'*[y_ex;u_ex] + (1/2)*([y_ex;u_ex]'*(Q*[y_ex;u_ex])) + obj_const_term
        obj_gap = solution_struct.obj_val - ex_obj_val;
        fprintf(['The objective value is %.2e. The gap w.r.t. the optimal objective is:\n',...
                 'Objective Gap = %.2e.\n'],solution_struct.obj_val, obj_gap);

        fprintf(['The relative errors w.r.t. the optimal solution are: \n', ...
                '|y-y_ex|/|y_ex| = %.2e, |u-u_ex|/|u_ex| = %.2e, |v-v_ex|/|v_ex| = %.2e.\n'],...
                y_relerr, u_relerr, p_relerr);         
        
        solution_statistics.y_ex = y_ex;
        solution_statistics.u_ex = u_ex;
        solution_statistics.p_ex = p_ex;
        solution_statistics.ex_obj_val = ex_obj_val;
    elseif (problem_choice == 3)
        % ============================================================================================== %
        % Convection Diffusion with H^1-regularization and state and control bounds.
        % ---------------------------------------------------------------------------------------------- %   
        problem_name = "Bounded_Convection_Diffusion";
        if (default == 1)
            epsilon = 0.01;                 % diffusion coefficient
            y_alpha = -0.1;  y_beta = 1;       % Upper and lower bounds for state and control.
            u_alpha = -0.6; u_beta = 0.6; 
        end   
        O = sparse(np,np);
        % Compute connectivity, stiffness, convection and mass matrices
        [ev,ebound] = q1grid(x_1x_2,mv,bound,mbound);
        [K,N,J_y,epe,eph,epw] = femq1_cd(x_1x_2,ev);

        % Compute SUPG stabilization matrix
        epe = epe/epsilon;
        esupg = find(epe <= 1); expe = epe;
        if any(expe)
           supg = inf;
           if isinf(supg)
              expe = 0.5*(1-1./expe);
              expe(esupg) = inf;
           else
              expe = ones(size(expe)); expe = supg*expe; expe(esupg) = inf;
           end
           epp = expe; epp(esupg) = 0; epp = epp.*eph./epw;
           S = femq1_cd_supg(x_1x_2,ev,expe,eph,epw);
        end

        % Compute relevant matrices for optimization algorithm
        D = epsilon*K+N+S; J_u = J_y+K;

        % Specify vectors relating to desired state and Dirichlet BCs
        yhat_vec = x_1x_2(:,1).^2.*x_1x_2(:,2).^2.*(x_1x_2(:,1) <= 0).*(x_1x_2(:,2) <= 0);
        bc_nodes = x_1x_2(bound,1).^2.*x_1x_2(bound,2).^2.*(x_1x_2(bound,1) <= 0).* ...
                   (x_1x_2(bound,2) <= 0);

        % Initialize RHS vector corresponding to desired state
        Jyhat = J_y*yhat_vec;

        % Enforce Dirichlet BCs on state, and zero Dirichlet BCs on adjoint
        [D,b] = nonzerobc_input(D,zeros(np,1),x_1x_2,bound,bc_nodes);
        [J_y,Jyhat] = nonzerobc_input(J_y,Jyhat,x_1x_2,bound,bc_nodes);
        J_constr = J_y; for i = 1:length(bound), J_constr(bound(i),bound(i)) = 0; end
        
        
        obj_const_term = (1/2).*(yhat_vec'*Jyhat);      % Constant term included in the objective.
        c = [-Jyhat; zeros(np,1)];
    
        
        A = [D J_constr];
        Q = [J_y O; O gamma.*J_u];
        lb = y_alpha.*ones(2*np,1);
        ub = y_beta.*ones(2*np,1);
        lb(np+1:2*np) = u_alpha.*ones(np,1);
        ub(np+1:2*np) = u_beta.*ones(np,1);
        [solution_struct] = Set_Up_IP_PMM(A,Q,b,c,obj_const_term,lb,ub,tol,IPM_maxit,printlevel,la_mode,fid);

        % Reasonable bound constraints are as follows:
        % epsilon = 1/100, beta = 1:      y \in [ 0   , 1], u \in [-0.007 , 0.007]
        % epsilon = 1/100, beta = 10^-1:  y \in [ 0   , 1], u \in [-0.17  , 0.12]
        % epsilon = 1/100, beta = 10^-2:  y \in [-0.1 , 1], u \in [-0.6   , 0.6]
        % epsilon = 1/100, beta = 10^-3:  y \in [-0.6 , 1], u \in [-0.6   , 1]
        % epsilon = 1/100, beta = 10^-4:  y \in [ 0   , 1], u \in [-0.8   , 1.2]
        % epsilon = 1/100, beta = 10^-5:  y \in [ 0   , 1], u \in [-1     , 1.2]
        % epsilon = 1/500, beta = 1:      y \in [ 0   , 1], u \in [-0.009 , 0.009]
        % epsilon = 1/500, beta = 10^-1:  y \in [-0.1 , 1], u \in [-0.35  , 0.3]
        % epsilon = 1/500, beta = 10^-2:  y \in [-0.8 , 1], u \in [-1.8   , 2]
        % epsilon = 1/500, beta = 10^-3:  y \in [-0.1 , 1], u \in [-1.2   , 1.4]
        % epsilon = 1/500, beta = 10^-4:  y \in [ 0   , 1], u \in [-0.9   , 1.2]
        % epsilon = 1/500, beta = 10^-5:  y \in [ 0   , 1], u \in [-1.2   , 1.2]
        % ______________________________________________________________________________________________ %

    elseif (problem_choice == 4)
        problem_name = "Bounded_Poisson_L1_Control";
        if (default == 1)
            u_alpha = -2; % control is constrained to be within range [-alpha,alpha]
            u_beta = 1.5;
            beta = 1e-6;    % regularization parameter of the L1 norm
        end

        O = sparse(np,np);

        % Compute connectivity, stiffness and mass matrices (D and J)
        [ev,ebound] = q1grid(x_1x_2,mv,bound,mbound);
        [D,J_y] = femq1_diff(x_1x_2,ev);
        R = (sum(J_y))'; % equivalently diag(sum(J)), as J is symmetric
        R(bound) = 0;  % account for the boundary conditions
        % Specify vectors relating to desired state, source term and Dirichlet BCs
        %yhat_vec = 1-(1+4*gamma*pi^4)/(2*pi^2)*sin(pi*x_1x_2(:,1)).*sin(pi*x_1x_2(:,2));

        yhat_vec = sin(pi*x_1x_2(:,1)).*sin(pi*x_1x_2(:,2));
        bc_nodes = ones(length(bound),1);

        % Initialize RHS vector corresponding to desired state
        Jyhat = J_y*yhat_vec;

        % Enforce Dirichlet BCs on state, and zero Dirichlet BCs on adjoint
        [D,b] = nonzerobc_input(D,zeros(np,1),x_1x_2,bound,bc_nodes);
        [J_y,Jyhat] = nonzerobc_input(J_y,Jyhat,x_1x_2,bound,bc_nodes);
        J_constr = J_y; for i = 1:length(bound), J_constr(bound(i),bound(i)) = 0; end
        obj_const_term = (1/2).*(yhat_vec'*Jyhat);      % Constant term included in the objective.
        c = [-Jyhat; (beta.*R).*ones(np,1); (beta.*R).*ones(np,1)];
    
        J_u = gamma.*J_y;
        A = [D -J_constr J_constr];
        Q = [J_y     O      O;
             O      J_u  -J_u;
             O     -J_u   J_u];
        diag_J_u = spdiags(J_u,0);
        diag_J_y = spdiags(J_y,0);
        Q_tilde = [spdiags(diag_J_y,0,np,np)     O                             O;
                   O        spdiags(diag_J_u,0,np,np)  -spdiags(diag_J_u,0,np,np);
                   O       -spdiags(diag_J_u,0,np,np)   spdiags(diag_J_u,0,np,np)];
        lb = -Inf.*ones(3*np,1);
        ub = Inf.*ones(3*np,1);
        lb(np+1:2*np) = max(u_alpha,0).*ones(np,1);
        ub(np+1:2*np) = max(u_beta,0).*ones(np,1);
        lb(2*np+1:3*np) = -min(u_beta,0).*ones(np,1);
        ub(2*np+1:3*np) = -min(u_alpha,0).*ones(np,1);
        [solution_struct] = Set_Up_IP_PMM(A,Q,b,c,obj_const_term,lb,ub,tol,IPM_maxit,printlevel,la_mode,fid,Q_tilde);

    elseif (problem_choice == 5)
        problem_name = "Conv_Diff_L1_Control";
        if (default == 1)
            u_alpha = -2; % control is constrained to be within range [-alpha,alpha]
            u_beta = 1.5;
            y_alpha = -inf;  y_beta = inf;       % Upper and lower bounds for state and control.
            beta = 1e-4;    % regularization parameter of the L1 norm
            epsilon = 0.02; % diffusion coefficient
        end

        O = sparse(np,np);

        % Compute connectivity, stiffness and mass matrices (D and J)
        [ev,ebound] = q1grid(x_1x_2,mv,bound,mbound);
        [K,N,J_y,epe,eph,epw] = femq1_cd(x_1x_2,ev);

        % Compute SUPG stabilization matrix
        epe = epe/epsilon;
        esupg = find(epe <= 1); expe = epe;
        if any(expe)
           supg = inf;
           if isinf(supg)
              expe = 0.5*(1-1./expe);
              expe(esupg) = inf;
           else
              expe = ones(size(expe)); expe = supg*expe; expe(esupg) = inf;
           end
           epp = expe; epp(esupg) = 0; epp = epp.*eph./epw;
           S = femq1_cd_supg(x_1x_2,ev,expe,eph,epw);
        end

        % Compute relevant matrices for optimization algorithm
        D = epsilon*K+N+S;
        
        R = (sum(J_y))'; % equivalently diag(sum(J)), as J is symmetric
        R(bound) = 0;  % account for the boundary conditions


        % Specify vectors relating to desired state and Dirichlet BCs
        yhat_vec = exp(-64.*((x_1x_2(:,1)-0.5).^2 + (x_1x_2(:,2)-0.5).^2));
        bc_nodes = 0.*x_1x_2(bound,1);
        
        
      %  yhat_vec = x_1x_2(:,1).^2.*x_1x_2(:,2).^2.*(x_1x_2(:,1) <= 0).*(x_1x_2(:,2) <= 0);
       % bc_nodes = x_1x_2(bound,1).^2.*x_1x_2(bound,2).^2.*(x_1x_2(bound,1) <= 0).* ...
       %            (x_1x_2(bound,2) <= 0);

        % Initialize RHS vector corresponding to desired state
        Jyhat = J_y*yhat_vec;

        % Enforce Dirichlet BCs on state, and zero Dirichlet BCs on adjoint
        [D,b] = nonzerobc_input(D,zeros(np,1),x_1x_2,bound,bc_nodes);
        [J_y,Jyhat] = nonzerobc_input(J_y,Jyhat,x_1x_2,bound,bc_nodes);
        J_constr = J_y; for i = 1:length(bound), J_constr(bound(i),bound(i)) = 0; end
        obj_const_term = (1/2).*(yhat_vec'*Jyhat);      % Constant term included in the objective.
        c = [-Jyhat; (beta.*R).*ones(np,1); (beta.*R).*ones(np,1)];
    
        J_u = gamma.*J_y;
        A = [D -J_constr J_constr];
        Q = [J_y     O      O;
             O      J_u  -J_u;
             O     -J_u   J_u];
        diag_J_u = spdiags(J_u,0);
        diag_J_y = spdiags(J_y,0);
        Q_tilde = [spdiags(diag_J_y,0,np,np)     O                             O;
                   O        spdiags(diag_J_u,0,np,np)  -spdiags(diag_J_u,0,np,np);
                   O       -spdiags(diag_J_u,0,np,np)   spdiags(diag_J_u,0,np,np)];
        lb = -Inf.*ones(3*np,1);
        ub = Inf.*ones(3*np,1);
        lb(1:np) = y_alpha.*ones(np,1);
        ub(1:np) = y_beta.*ones(np,1);
        lb(np+1:2*np) = max(u_alpha,0).*ones(np,1);
        ub(np+1:2*np) = max(u_beta,0).*ones(np,1);
        lb(2*np+1:3*np) = -min(u_beta,0).*ones(np,1);
        ub(2*np+1:3*np) = -min(u_alpha,0).*ones(np,1);
        [solution_struct] = Set_Up_IP_PMM(A,Q,b,c,obj_const_term,lb,ub,tol,IPM_maxit,printlevel,la_mode,fid,Q_tilde);

    end
    solution_statistics.y = solution_struct.x(1:np);
    if (problem_choice <= 3)
        solution_statistics.u = solution_struct.x(np+1:2*np,1);
    else
        solution_statistics.u = solution_struct.x(np+1:2*np,1) - solution_struct.x(2*np+1:3*np,1);
    end
    solution_statistics.v = solution_struct.y;
    solution_statistics.z_l = solution_struct.z_l;
    solution_statistics.z_u = solution_struct.z_u;
    solution_statistics.total_time = solution_struct.pre_time + solution_struct.runtime;
    solution_statistics.objective_value = solution_struct.obj_val;
    solution_statistics.status = solution_struct.opt;
    solution_statistics.total_Krylov_iters = solution_struct.Krylov_its;
    solution_statistics.max_nnzL = solution_struct.max_nnzL;
    solution_statistics.total_iters = solution_struct.IP_iter;
    if (solution_struct.opt == 1)                                       % Success
       solution_statistics.problem_converged = 1;
       fprintf(fileID,'The optimal solution objective is %d.\n',solution_struct.obj_val);
    elseif (solution_struct.opt == 2)                                   % Primal Infeasibility
        fprintf(fileID,['Primal infeasibility has been detected.\n',...
                        'Returning the last iterate.\n']); 
    elseif (solution_struct.opt == 3)                                   % Dual Infeasibility
        fprintf(fileID,['Dual infeasibility has been detected\n',...
                        'Returning the last iterate.\n']); 
    elseif (solution_struct.opt == 4)                                   % Numerical Error
        fprintf(fileID,['Method terminated due to numerical error.\n',...   
                        'Returning the last iterate.\n']); 
    else                                                                % Reached maximum iterations
       fprintf(fileID,['Maximum number of iterations reached.\n',...
                           'Returning the last iterate.\n']); 
    end
    fprintf(fileID,['Name = %s & IP iters = %d  & Krylov iters = %d & max_nnzL = %.2e & ' ...
        'Time = %.2e & opt = %s  \n'],problem_name, solution_struct.IP_iter, solution_struct.Krylov_its,...
                                     solution_struct.max_nnzL, solution_statistics.total_time, string(solution_struct.opt == 1)); 
end

