function [solution_struct] = Set_Up_IP_PMM(A,Q,b,c,obj_const_term,lb,ub,tol,maxit,printlevel,la_mode,fid,Q_tilde)
% ==================================================================================================================== %
% This function is prepares the required structs to call an Interior Point-Proximal Method of Multipliers, 
% suitable for solving linear and convex quadratic
% programming problems. The method takes as input a problem of the following form:
%
%                                    min   c^T x + (1/2)x^TQx,
%                                    s.t.  A x = b,
%                                          lb <= x <= ub    x in R^n, A in R^(m)x(n),
%
% and solves it to optimality, returning the primal and dual optimal solutions (or a message indicating that the
% optimal solution was not found). The solution statistics are gathered in a struct and returned by this function.
%
% INPUT PARAMETERS/ALLOWED INPUT FORMATS:
%
%    Set_Up_IP_PMM(A, Q, b, c): Basic problem data 
%                               A -> mxn constraint matrix (if empty use sparse mxn matrix),
%                               Q -> nxn Hessian matrix (if empty use sparse nxn matrix),
%                               b -> the right hand side vector (if empty use mx1 zero vector),
%                               c -> the linear part of the objective (if empty use nx1 zero vector).
%                               By default we then set all remaining parameters to empty
%                               and they are initialized within IP-PMM (see the relevant documentation).
%    Set_Up_IP_PMM(A, Q, b, c, obj_const_term): An input of any constant term missing from the objective due to  
%                                               modelling choices.
%    Set_Up_IP_PMM(A, Q, b, c, obj_const_term, lb, ub): The user may add lower and upper bounds respectively in that way.
%    Set_Up_IP_PMM(A, Q, b, c, obj_const_term, lb, ub, tol): The user may set the required tolerance in that way.
%    Set_Up_IP_PMM(A, Q, b, c, obj_const_term, lb, ub, tol, maxit): Adds maximum number of IP-PMM iterations.
%    Set_Up_IP_PMM(A, Q, b, c, obj_const_term, lb, ub, tol, maxit, printlevel): sets the printlevel.
%                                                           0: turn off iteration output
%                                                           1: print primal and dual residual and duality measure
%                                                           2: print step length
%                                                           3: print intermediate Krylov iterates
%    Set_Up_IP_PMM(A, Q, b, c, obj_const_term, lb, ub, tol, maxit, printlevel, la_mode):
%                                                           la_mode (linear algebra mode):
%                                                           "exact" -> factorization based code (default),
%                                                           "inexact" -> iterative-based code.
%           "inexact" is usually better than "exact" on larger problems, and it requires less memory.
%    Set_Up_IP_PMM(A,Q,b,c,obj_const_term,lb,ub,tol,maxit,printlevel,la_mode,fid): fid sets 
%                                                                                  the printing file.
%
% OUTPUT STRUCT:
% solution_struct with fields: x (primal solution), y (dual multipliers), z_l (dual slacks corresponding to 
%                              the lower bounded variables), z_u (dual slacks corr. to upper bounded var.),
%                              obj_val (the value of the objective function at the terminating iterate of IP-PMM).
%
% Author: Spyridon Pougkakiotis, January 2021, Edinburgh.
% ____________________________________________________________________________________________________________________ %
    if (nargin < 11 || isempty(la_mode))
        la_mode = "exact";
    end
    if (nargin < 12 || isempty(fid))
        fid = 1;
    end
    tic;
    % ================================================================================================================ %
    % Convert the problem to the IP-PMM format. We have two choices (the former being preferable):
    %           1) Non-negative and free variables: Each variable x_i is either free, or x_i >= 0.
    %           2) each variable can only have either a lower or an upper bound (abritrary bounds).
    % ---------------------------------------------------------------------------------------------------------------- %
    n_prev = size(c,1);
    m_prev = size(b,1);
    format_option = "semi_standard_form";
    if (format_option == "semi_standard_form")
        [A, Q, b, c, lb, ub, obj_const_term_2] = QP_Convert_to_Semi_Standard_Form(A, Q, b, c, lb, ub);
    elseif (format_option == "lower_or_upper_bound_form")
        [A, Q, b, c, lb, ub, obj_const_term_2] = QP_Convert_to_L_Or_U_Form(A, Q, b, c, lb, ub);
    end
    m = size(b,1);
    n = size(c,1);  
    if (nargin < 13 || isempty(Q_tilde))
        Q_tilde = Q;
    elseif (n_prev < n)
        Q_tilde = [Q_tilde                   sparse(n_prev,n-n_prev);
                   sparse(n-n_prev,n_prev)   sparse(n-n_prev,n-n_prev)];
    end
    % ________________________________________________________________________________________________________________ %
    
    % ================================================================================================================ %
    % Scale the problem matrix
    % ---------------------------------------------------------------------------------------------------------------- %
    scaling_direction = 'l';
    scaling_mode = 0;
     if (scaling_direction == 'r') % Apply the right scaling.
        [D] = Scale_the_problem(A,scaling_mode,scaling_direction);
        A = A*spdiags(D,0,n,n); 
        c = c.*D;
        Q = spdiags(D,0,n,n)*Q*spdiags(D,0,n,n);
        lb = lb./D;
        ub = ub./D;
    elseif (scaling_direction == 'l') % Apply the left scaling.
        [D] = Scale_the_problem(A,scaling_mode,scaling_direction);
        A = spdiags(D,0,m,m)*A;                         
        b = b.*D;
     elseif (scaling_direction == 'b') % Apply left and right scaling
        [D_R] = Scale_the_problem(A,scaling_mode,'r');
        [D_L] = Scale_the_problem(A,scaling_mode,'l');
        A = spdiags(D_L,0,m,m)*A*spdiags(D_R,0,n,n);
        b = b.*D_L;
        c = c.*D_R;
        Q = spdiags(D_R,0,n,n)*Q*spdiags(D_R,0,n,n);
        lb = lb./D_R;
        ub = ub./D_R;
    end
    % ________________________________________________________________________________________________________________ %
    lb_vars = (lb > -Inf);                                    % logic 1 for components of x having a lower bound.
    ub_vars = (ub < Inf);                                     % logic 1 for components of x having an upper boud.
    A_tr = A';
    solution_struct = struct();
    A_struct = struct();
    Q_struct = struct();
    Newton_struct = struct();
    if (la_mode == "inexact") Precond_struct = struct(); end
    % ================================================================================================================ %
    % Function handles and measures concerning matrix A.
    % ---------------------------------------------------------------------------------------------------------------- %
    A_struct.A = @(w) A*w;
    A_struct.A_tr = @(w) A_tr*w;
    A_struct.NE_diag = sum(A.^2,2);
    A_struct.A_1_norm = norm(A,1);
    A_struct.A_inf_norm = norm(A,'inf');
    % ________________________________________________________________________________________________________________ %
    
    % ================================================================================================================ %
    % Function handles and measures concerning matrix Q.
    % ---------------------------------------------------------------------------------------------------------------- %
    Q_struct.Q = @(w) Q*w;
    Q_struct.diag = spdiags(Q,0);
    Q_struct.nnz = nnz(Q);
    Q_struct.Q_1_norm = norm(Q,1);
    Q_struct.Q_inf_norm = norm(Q,'inf');
    % ________________________________________________________________________________________________________________ %

    % ================================================================================================================ %
    % Function handles and information required to build and solve the Newton system.
    % ---------------------------------------------------------------------------------------------------------------- %
    Newton_struct.la_mode = la_mode;
    Newton_struct.IR = "inner";                        % If IR == "outer" -> iterative refinement in Newton solution,
                                                       % If IR == "inner" -> iterative refinement in preconditioner,
                                                       % If IR == "inner_outer" -> both IR methods at once,
                                                       % If IR == "none"  -> no iterative refinement.
    if (la_mode == "exact")
        % Set the factorization function.
        Newton_struct.Newton_fact = @(x,z_l,z_u,delta,rho,pivot_thr) ...
                                  Newton_factorization(A,A_tr,Q,x,z_l,z_u,delta,rho,lb,ub,lb_vars,ub_vars,pivot_thr,Newton_struct.IR);
        % Set the backsolve function.
        Newton_struct.Newton_backsolve = @(NSdata,res_p,res_d,res_mu_l,res_mu_u) ...
                                      Newton_backsolve(NSdata,res_p,res_d,res_mu_l,res_mu_u);
    elseif (la_mode == "inexact")
        Newton_struct.prec_approach = "Chol_based_preconditioner";   % Use "Chol_based_preconditioner" for Schur complement approximation,
                                                                     % or "LDL_based_preconditioner" for LDL'-based approximation.
        if (Newton_struct.prec_approach == "Chol_based_preconditioner") 
            double_approximation = true;                       % If we use Cholesky-based preconditioner, drop dense rows and columns.  
        else
            double_approximation = false;
            Newton_struct.minres_setting = "Schur_complement";     % Use "augmented_system" for L|D|L^T preconditioner,
                                                                    % or "Schur_complement" for implicit Schur complement
                                                                    % approximation.
        end
        % ============================================================================================================ %
        % Drop possible dense rows within the preconditioner (only in the inexact case).
        % ------------------------------------------------------------------------------------------------------------ %
        k_max = 15;
        dens_c = 0.15;
        dens_r = 0.4;
        if (double_approximation && k_max > 0)
            [R, C] = find_row_and_col_drop(k_max, A, A_tr, m, n, dens_c, dens_r);      
        else
            R = true(m,1);
            C = true(n,1);
        end
        % ____________________________________________________________________________________________________________ %
        % ============================================================================================================ %
        % Set the Newton system data, operators (normal equations + augmented system), and the Krylov solver.
        % ------------------------------------------------------------------------------------------------------------ %
        Newton_struct.AS_multiplier = @ (x,NSdata) AS_multiplier(x,NSdata);
        Newton_struct.NE_multiplier = @(x,NSdata)  NE_multiplier(x,NSdata);
        Newton_struct.Newton_itersolve = @ (pred,NSdata,PS,res_p,res_d,res_mu_l,res_mu_u,maxit,tol) ...
                            Newton_itersolve(fid,pred,NSdata,PS,res_p,res_d,res_mu_l,res_mu_u,maxit,tol,printlevel);
        Newton_struct.Newton_matrix_data = @ (iter,mu,A_struct,Q_struct,x,z_l,z_u,delta,rho,pivot_thr,solver,prec_approach,p_inf,d_inf) ...
                            Newton_matrix_data(iter,mu,A_struct,Q_struct,x,z_l,z_u,lb,ub,lb_vars,ub_vars,delta,rho,m,n,...
                                               pivot_thr,solver,prec_approach,Newton_struct.IR,p_inf,d_inf);
        % ____________________________________________________________________________________________________________ %
        % ============================================================================================================ %
        % Set the preconditioner construction, tuning, and backolsve functions.
        % ------------------------------------------------------------------------------------------------------------ %        
        Precond_struct.build_out_precond = @ (NSdata,droptol,mu,set_proximal_inner_IR) ...
                            build_outer_preconditioner(A,A_tr,Q_tilde,NSdata,droptol,mu,R,C,double_approximation,set_proximal_inner_IR);
        Precond_struct.set_out_precond_tol = @ (droptol,iter,itmax,nnzL,roof) ...
                            precond_settol(droptol,iter,itmax,nnzL,roof);
        Precond_struct.out_precond_backsolve = @ (x,PS,solver) Precond_Operator(x,PS,solver);
        Newton_struct.Precond_struct = Precond_struct;
        % ____________________________________________________________________________________________________________ %
    end
    % ________________________________________________________________________________________________________________ %
    solution_struct.pre_time = toc;
    tic;
    [x,y,z_l,z_u,opt,IP_iter,Krylov_its,max_nnzL] = IP_PMM(A_struct,Q_struct,Newton_struct,b,c,lb,ub,tol,maxit,printlevel,fid);
    solution_struct.runtime = toc;
    solution_struct.max_nnzL = max_nnzL;
    solution_struct.x = x; solution_struct.y = y; solution_struct.z_l = z_l;
    solution_struct.z_u = z_u; solution_struct.IP_iter = IP_iter; solution_struct.opt = opt;
    solution_struct.obj_val = c'*x + (1/2)*(x'*(Q*x)) + obj_const_term + obj_const_term_2;
    solution_struct.Krylov_its = Krylov_its;
end

