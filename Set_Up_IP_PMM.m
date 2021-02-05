function [solution_struct] = Set_Up_IP_PMM(A,Q,b,c,obj_const_term,lb,ub,tol,maxit,printlevel)
% ==================================================================================================================== %
% This function is prepares the required structs to call an Interior Point-Proximal Method of Multipliers, 
% suitable for solving linear and convex quadratic
% programming problems. The method takes as input a problem of the following form:
%
%                                    min   c^T x + (1/2)x^TQx,
%                                    s.t.  A x = b,
%                                          lb <= x <= ub,           x in R^n, A in R^(m)x(n),
%
% and solves it to optimality, returning the primal and dual optimal solutions (or a message indicating that the
% optimal solution was not found). The solution statistics are gathered in a struct and returned by this function.
%
% INPUT PARAMETERS/ALLOWED INPUT FORMATS:
%
%    Set_Up_IP_PMM(A, Q, b, c): Basic problem data 
%                               A -> mxn constraint matrix (can be empty),
%                               Q -> nxn Hessian matrix (can be empty),
%                               b -> the right hand side vector (cannot be empty),
%                               c -> the linear part of the objective (cannot be empty).
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
% OUTPUT STRUCT:
% solution_struct with fields: x (primal solution), y (dual multipliers), z_l (dual slacks corresponding to 
%                              the lower bounded variables), z_u (dual slacks corr. to upper bounded var.),
%                              obj_val (the value of the objective function at the terminating iterate of IP-PMM).
%
% Author: Spyridon Pougkakiotis, January 2021, Edinburgh.
% ____________________________________________________________________________________________________________________ %
    A_tr = A';
    m = size(b,1); n = size(c,1);
    solution_struct = struct();
    A_struct = struct();
    Q_struct = struct();
    Newton_struct = struct();
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
    Q_struct.diag = spdiags(Q,0,n,n);
    Q_struct.Q_1_norm = norm(Q,1);
    Q_struct.Q_inf_norm = norm(Q,'inf');
    % ________________________________________________________________________________________________________________ %
    
    % ================================================================================================================ %
    % Function handles and information required to build and solve the Newton system.
    % ---------------------------------------------------------------------------------------------------------------- %
    Newton_struct.Newton_fact = @(x,z_l,z_u,delta,rho,lb,ub,lb_vars,ub_vars,pivot_thr) ...
                                  Newton_factorization(A,A_tr,Q,x,z_l,z_u,delta,rho,lb,ub,lb_vars,ub_vars,pivot_thr);
    Newton_struct.Newton_backsolve = @(NS,res_p,res_d,res_mu_l,res_mu_u) ...
                                       Newton_backsolve(NS,res_p,res_d,res_mu_l,res_mu_u);
    % ________________________________________________________________________________________________________________ %

    [x,y,z_l,z_u,opt,IP_iter] = IP_PMM(A_struct,Q_struct,Newton_struct,b,c,lb,ub,tol,maxit,printlevel);
    solution_struct.x = x; solution_struct.y = y; solution_struct.z_l = z_l;
    solution_struct.z_u = z_u; solution_struct.IP_iter = IP_iter; solution_struct.opt = opt;
    solution_struct.obj_val = c'*x + (1/2)*(x'*(Q*x)) + obj_const_term;
end

