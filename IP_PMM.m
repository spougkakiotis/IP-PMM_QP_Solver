function [x,y,z_l,z_u,opt,IP_iter,Krylov_its,max_nnzL] = IP_PMM(A_struct,Q_struct,Newton_struct,b,c,lb,ub,tol,maxit,printlevel,fid)
% ==================================================================================================================== %
% This function is an Interior Point-Proximal Method of Multipliers, suitable for solving linear and convex quadratic
% programming problems. The method takes as input a problem of the following form:
%
%                                    min   c^T x + (1/2)x^TQx,
%                                    s.t.  A x = b,
%                                          lb <= x <= ub,           x in R^n, A in R^(m)x(n),
%
% and solves it to optimality, returning the primal and dual optimal solutions (or a message indicating that the
% optimal solution was not found).
%
% INPUT PARAMETERS/ALLOWED INPUT FORMATS:
% IP_PMM(A, Q, Newton_struct, b, c): Basic problem data 
%                     A_struct -> Constraint matrix 
%                                 must be given as a struct of function handles, including:
%                                 A_struct.A(w) = A*w,
%                                 A_struct.A_tr(w) = A'*w,
%                                 A_struct.NE_diag(D) = diag(A*Diag(D)*A'), where D is a vector.
%                                 If empty, the default operations become the zero function (sets all inputs to zero).
%                     Q_struct -> Hessian matrix
%                                 must be given as a struct of function handles, including: 
%                                 Q_struct.Q(w) = Q*w,
%                                 Q_struct.diag = diag(Q).
%                                 If empty, the default operations become the zero function (sets all inputs to zero).
%                     Newton_struct -> Contains all relevant function handles needed to solve the Newton system.
%                                      It cannot be empty. For more on the functions, see their respective 
%                                      documentation.
%                     b -> the right hand side of the equality constraints.
%                          **It cannot be empty**. 
%                     c -> the linear coefficients of the objective function.
%                          **It cannot be empty**. 
%
%                     ------------------------------------------------------------------------
%                     Default lb = -Inf.*ones(n,1),  Default ub = Inf.*ones(n,1)
%                     (**** Without bound constraints the method becomes a standard PMM. ****)
%                     ________________________________________________________________________
%                     Default tolerance = 1e-6, Default maximum IP-PMM iterations = 100.
%
% IP_PMM(A, Q, Newton_struct, b, c, lb, ub): This way, the user can specify potential bound constraints on x
%                             lb -> vector of size n, containing lower bounds on x 
%                                   Set lb(j) = -Inf, if x(j) does not have a lower bound.
%                             ub -> vector of size n, containing the upper bounds on x
%                                   Set ub(j) = Inf, if x(j) does not have an upper bound.
%
% IP_PMM(A, Q, Newton_struct, b, c, lb, ub, tol): Specify the tolerance to which the problem is solved.
%
% IP_PMM(A, Q, Newton_struct, b, c, lb, ub, tol, max_it): Specify the maximum number of IP-PMM iterations.
%
% IP_PMM(A, Q, Newton_struct, b, c, lb, ub, tol, max_it, printlevel): sets the printlevel.
%                                              0: turn off iteration output
%                                              1: print primal and dual residual and duality measure
%                                              2: print step length
%                                              3: print intermediate Krylov iterate statistics
% IP_PMM(A, Q, Newton_struct, b, c, lb, ub, tol, max_it, printlevel,fid): prints on a file with indicator fid.
% OUTPUT: [x,y,z,opt,iter], where:
%         x: primal solution
%         y: Lagrange multiplier vector
%         z_l: dual slack variables corresponding to the lower bound constraints
%         z_u: dual slack variables corresponding to the upper bound constraints
%         opt: the possible outcomes are the following:
%                   0: non-optimal -> maximum iterations reached
%                   1: optimal -> converge to the accepted tolerance
%                   2: infeasible -> primal infeasible 
%                   3: infeasible -> dual infeasible
%                   4: ill-conditioning -> problem terminated due to numerical error
%                   5: inaccuracy -> problem terminated due to insufficient accuracy
%         iter: numeber of iterations to termination.
%
% Author: Spyridon Pougkakiotis, January 2021, Edinburgh.
% ____________________________________________________________________________________________________________________ %

% ==================================================================================================================== %
% Parameter filling and dimensionality testing.
% -------------------------------------------------------------------------------------------------------------------- %
if (isempty(b) || isempty(c))
    error("The right-hand side, b, and the linear objective coefficients, c, cannot be empty.\n");
end
if (issparse(b))  b = full(b);   end
if (issparse(c))  c = full(c);   end

% Make sure that b and c are column vectors of dimension m and n.
if (size(b,2) > 1) b = (b)'; end
if (size(c,2) > 1) c = (c)'; end
m = size(b,1); n = size(c,1);

zero_function = @(x) zeros(size(x,1),1);
if (isempty(A_struct)) 
    A_struct.A = zero_function;
    A_struct.A_tr = zero_fuction;
    A_struct.NE_diag = zeros(m,1);
end
if (isempty(Q_struct))
    Q_struct.Q = zero_function; 
    Q_struct.diag = zeros(n,1);
end

% Make sure that A and Q are given as function handles.
if (~isa(A_struct,'struct')) 
    error(["The constraint matrix A should be given as a struct of function handles.\n",...
           "For more information, see the description of the algorithm."]); 
end
if (~isa(Q_struct,'struct')) 
    error(["The Hessian matrix Q should be given as a function handle.\n",...
           "For more information, see the description of the algorithm."]); 
end
% Set default values for missing parameters.
if (nargin < 6  || isempty(lb))            lb = -Inf.*ones(n,1);          end
if (nargin < 7  || isempty(ub))            ub = Inf.*ones(n,1);           end
if (nargin < 8 || isempty(tol))            tol = 1e-4;                    end
if (nargin < 9 || isempty(maxit))          maxit = 100;                   end
if (nargin < 10 || isempty(printlevel))    printlevel = 1;                end
if (nargin < 11 || isempty(fid))           fid = 1;                       end
% ____________________________________________________________________________________________________________________ %

% ==================================================================================================================== %
% Set indexes for acting bound constraints. 
% -------------------------------------------------------------------------------------------------------------------- %
lb_vars = (lb > -Inf);                                    % logic 1 for components of x having a lower bound.
ub_vars = (ub < Inf);                                     % logic 1 for components of x having an upper boud.
num_lb_vars = nnz(lb_vars); num_ub_vars = nnz(ub_vars);   % number of variables with lower (upper bound, resp.)
% ____________________________________________________________________________________________________________________ %

% ==================================================================================================================== %
% Starting point (Set mu = 0 if there are no bound constraints at all).
% -------------------------------------------------------------------------------------------------------------------- %
[x,y,z_l,z_u] = Mehrotra_Warm_Start(A_struct, Q_struct, b, c, lb, ub, lb_vars, ub_vars);

lambda = y;   zeta = x;                                   % Initial estimate of the primal and dual optimal solution.
if (num_lb_vars + num_ub_vars > 0)                        % Initial value of mu (when bound constraints are present).
    mu = (x(lb_vars)-lb(lb_vars))'*z_l(lb_vars) + ...
         (ub(ub_vars) - x(ub_vars))'*z_u(ub_vars);
   % mu = max((mu)/(num_lb_vars + num_ub_vars),0.1);
    mu = (mu)/(num_lb_vars + num_ub_vars);
    res_mu_l = zeros(n,1); res_mu_u = zeros(n,1);         % To be used as residuals for complementarity.
    mu_prev = mu;                                         % Keep track the previous value of mu.
else                                                      % Switch to a pure PMM method (no inequality constraints).
    mu = 0;     
end
% ____________________________________________________________________________________________________________________ %

% ==================================================================================================================== %
% Parameter Initialization.
% -------------------------------------------------------------------------------------------------------------------- %
ctr_struct = build_control_structure();                   % Contains all the monitoring parameters of IP-PMM.
ctr_struct.isQP = (Q_struct.Q_1_norm > 0);                % isQP monitors whether we are solving a QP (true) or 
                                                          % an LP (false).
delta = 8;  rho = 8;                                      % Initial dual and primal regularization value.
reg_limit = max(0.1*min(tol,1e-6)*(1/max(1e0,max(A_struct.A_1_norm*A_struct.A_inf_norm,...
                Q_struct.Q_1_norm*Q_struct.Q_inf_norm))),5e-10);
reg_limit = min(reg_limit,1e-8);                          % Controlled perturbation.
pivot_thr = 1e-10;                                        % Minimum pivot threshold for factorization methods.  
if (Newton_struct.la_mode == "inexact")
    iterlin = 50;                                         % Estimate of Krylov iterations.
    droptol = 1e2;                                        % Initial dropping tolerance for the outer preconditioner.
    roof = 3e7;                                           % Dangerously large nnz of the preconditioner ->
                                                          % be more conservative.
    itertot = zeros(2*maxit,1);                           % Keep track the total number of inner and outer iterations.
    Krylov_tol = tol;                                     % Absolute tolerance for the Krylov solver.
    Krylov_its = 0;                                       % Keep track of total # of Krylov iterations.
    max_nnzL = 0;
    if (Q_struct.nnz == nnz(Q_struct.diag))  % Linear programming or Separable quadratic objective.
        solver = "pcg";
        maxit_Krylov = 100;
        fprintf(fid,"Employing pcg as Q is zero or diagonal. Normal equations are solved.\n"); 
    else
        solver = "minres";
        maxit_Krylov = 200;
        fprintf(fid,"Employing minres as Q is neither zero nor diagonal, or we prefer solving the augmented system.\n"); 
    end
end
IP_PMM_header(fid,printlevel,Newton_struct.la_mode);      % Set the header for printing.
% ____________________________________________________________________________________________________________________ %

while (ctr_struct.IP_iter < maxit)
% -------------------------------------------------------------------------------------------------------------------- %
% IP-PMM Main Loop structure:
% Until (||Ax_k - b|| < tol && ||c + Qx_k - A^Ty_k - z_k|| < tol && mu < tol) do
%   Choose sigma in [sigma_min, sigma_max] and solve:
%
%      [ -(Q + Theta^{-1} + rho I)   A^T            I]  (Delta x)    (c + Qx_k - A^T y_k - [z_C_k; 0] + rho (x-zeta))
%      [           A               delta I          0]  (Delta y)    (b - Ax_k - delta (y-lambda))
%      [         Z_C                 0            X_C]            =
%                                                       (Delta z)    (sigma e_C - X_C Z_C e_C)
%      [           0                 0            X_F]               (           0           )
%
%   where mu = x_C^Tz_C/|C|. Set (z_F)_i = 0, for all i in F.
%   Find two step-lengths a_x, a_z in (0,1] and update:
%              x_{k+1} = x_k + a_x Delta x, y_{k+1} = y_k + a_z Delta y, z_{k+1} = z_k + a_z Delta z
%   k = k + 1
% End
% -------------------------------------------------------------------------------------------------------------------- %
    if (ctr_struct.IP_iter > 1)
        nr_res_p = new_nr_res_p;
        nr_res_d = new_nr_res_d;
    else
        nr_res_p = b-A_struct.A(x);                                 % Non-regularized primal residual
        nr_res_d = c-A_struct.A_tr(y)-z_l + z_u + Q_struct.Q(x);    % Non-regularized dual residual.
    end
    res_p = nr_res_p - delta.*(y-lambda);                           % Regularized primal residual.
    res_d = nr_res_d + rho.*(x-zeta);                               % Regularized dual residual.
    p_inf = norm(nr_res_p,inf);
    d_inf = norm(nr_res_d,inf);
    % ================================================================================================================ %
    % Check terminatio criteria
    % ---------------------------------------------------------------------------------------------------------------- %
    if (p_inf/(1+norm(b,inf)) < tol && d_inf/(1+norm(c,inf)) < tol &&  mu < tol )
        fprintf('optimal solution found\n');
        ctr_struct.opt = 1;
        break;
    end
    if ((norm(y-lambda)> 10^10 && norm(res_p) < tol && ctr_struct.no_dual_update > 5)) 
        fprintf('The primal-dual problem is infeasible\n');
        ctr_struct.opt = 2;
        break;
    end
    if ((norm(x-zeta)> 10^10 && norm(res_d) < tol && ctr_struct.no_primal_update > 5))
        fprintf('The primal-dual problem is infeasible\n');
        ctr_struct.opt = 3;
        break;
    end
    % ________________________________________________________________________________________________________________ %
    ctr_struct.IP_iter = ctr_struct.IP_iter+1;
    if (num_lb_vars + num_ub_vars > 0)
        [x,z_l,z_u,mu] = boundary_control(x,z_l,z_u,lb,ub,lb_vars,ub_vars,num_lb_vars,num_ub_vars,mu);
        mu_prev = mu;
    end
    % ================================================================================================================ %
    % Avoid the possibility of converging to a local minimum -> Decrease the minimum regularization value.
    % ---------------------------------------------------------------------------------------------------------------- %
    if (reg_limit ~= 5e-13)
        [ctr_struct,reg_limit] = avoid_local_min(ctr_struct,rho,delta,reg_limit);
       % pivot_thr = reg_limit/5;
    end
    % ________________________________________________________________________________________________________________ %
    % ================================================================================================================ %
    % Compute the Newton factorization (exact case) or gather information about the Newton matrix (inexact case)
    % ---------------------------------------------------------------------------------------------------------------- %
    if (Newton_struct.la_mode == "exact")
        NSdata = Newton_struct.Newton_fact(x,z_l,z_u,delta,rho,pivot_thr);
    elseif (Newton_struct.la_mode == "inexact")
        NSdata = Newton_struct.Newton_matrix_data(ctr_struct.IP_iter,A_struct,Q_struct,x,z_l,z_u,...
                                                  delta,rho,pivot_thr,solver,Newton_struct.prec_approach,p_inf,d_inf); 
        if (Newton_struct.prec_approach == "LDL_based_preconditioner")   % Indicate which minres solver to use.
            NSdata.minres_setting = Newton_struct.minres_setting;
        end
    end
    % ________________________________________________________________________________________________________________ %
 
    % ================================================================================================================ %
    % Predictor step: Solve the Newton system and compute a centrality measure.
    % ---------------------------------------------------------------------------------------------------------------- %
    predictor = true;
    res_mu_l(lb_vars) = (lb(lb_vars)-x(lb_vars)).*z_l(lb_vars);
    res_mu_u(ub_vars) = (x(ub_vars)-ub(ub_vars)).*z_u(ub_vars);
    % ================================================================================================================ %
    % Solve the Newton system with the predictor right hand side -> Optimistic view, solve as if you wanted to 
    %                                                               solve the original problem in 1 iteration.
    % ---------------------------------------------------------------------------------------------------------------- %
    if (Newton_struct.la_mode == "exact")
        [dx,dy,dz_l,dz_u,instability] = Newton_struct.Newton_backsolve(NSdata,res_p,res_d,res_mu_l,res_mu_u); 
        nnzL = nnz(NSdata.L);
        if (ctr_struct.IP_iter == 1)
            max_nnzL = nnzL;
        end
    elseif (Newton_struct.la_mode == "inexact")
        if (ctr_struct.IP_iter > 1) 
            if (~PS.instability)
                if (~PS.double_approximation || PS.no_prec || Newton_struct.prec_approach == "LDL_based_preconditioner")
                    nnzL = nnz(PS.L_M);
                else
                    nnzL = nnz(PS.R_NE_PS.L_11) + nnz(PS.R_NE_PS.L_22);
                end
                if (nnzL > max_nnzL)
                    max_nnzL = nnzL;
                end
            end
        else
            nnzL = 0;
        end  
        [droptol] = Newton_struct.Precond_struct.set_out_precond_tol(droptol,iterlin,maxit_Krylov,nnzL,roof);
        PS = Newton_struct.Precond_struct.build_out_precond(NSdata,droptol,mu,ctr_struct.set_proximal_inner_IR); % Preconditioner for the Newton system.
        if (PS.instability == false)
            [dx,dy,dz_l,dz_u,instability,iterlin,drop_direction,Reinforce_Inner_IR] = ...
             Newton_struct.Newton_itersolve(predictor,NSdata,PS,res_p,res_d,res_mu_l,...
                                            res_mu_u,maxit_Krylov,Krylov_tol);
        else
            instability = true;
        end
        % Keep track of the number of iterations performed. If an iter. is repeated, also keep track of those iterations.
        itertot(2*ctr_struct.IP_iter-1) = itertot(2*ctr_struct.IP_iter-1) + iterlin;  
        if (drop_direction == true)
            [ctr_struct,mu,delta,rho] = control_direction_accuracy(ctr_struct,mu,mu_prev,delta,rho,predictor);
            if (ctr_struct.opt == 5)    % 5 signifies failure due to inaccuracy
                break;
            end
            continue;
        end
    end
    if (instability == true)        % Checking if the matrix is too ill-conditioned. Mitigate it.
        [ctr_struct,delta,rho,reg_limit] = control_ill_conditioning(ctr_struct,rho,delta,reg_limit,tol,predictor);
        if (ctr_struct.opt == 4)    % 4 signifies failure due to numerical instability.
            break;
        end
        continue;
    end
    ctr_struct.retry_p = 0;
    % ________________________________________________________________________________________________________________ %
      
    % ================================================================================================================ %
    % Step in the lower-or-upper bound (primal) and non-negativity orthant (dual).
    % ---------------------------------------------------------------------------------------------------------------- %
    if (num_lb_vars + num_ub_vars > 0)
        [alpha_x,alpha_z] = step_to_the_boundary(x,z_l,z_u,dx,dz_l,dz_u,lb,ub,lb_vars,ub_vars);
    else
        alpha_x = 1;         % If we have no inequality constraints, Newton method is exact -> Take full step.
        alpha_z = 1;        
    end    
    % ________________________________________________________________________________________________________________ %
    if (num_lb_vars + num_ub_vars > 0)   % Corrector only in the case we have inequality constraints.
        % ============================================================================================================ %
        % Corrector step: Solve Newton system with the corrector right hand side. Solve as if you wanted to direct the 
        %                 method in the center of the central path.
        % ------------------------------------------------------------------------------------------------------------ %
        predictor = false;
        centr_measure = (x(lb_vars)-lb(lb_vars) + alpha_x.*dx(lb_vars))'*(z_l(lb_vars) + alpha_z.*dz_l(lb_vars)) + ...
                             (ub(ub_vars) - x(ub_vars)-alpha_x.*dx(ub_vars))'*(z_u(ub_vars) + alpha_z.*dz_u(ub_vars));
        %mu = max((centr_measure/((num_lb_vars + num_ub_vars)*mu))^3*mu,0.1*mu_prev); 
        mu = (centr_measure/((num_lb_vars + num_ub_vars)*mu))^3*mu; 
      %  res_mu_l(lb_vars) = mu.*ones(num_lb_vars,1)-dx(lb_vars).*dz_l(lb_vars)-(x(lb_vars)-lb(lb_vars)).*z_l(lb_vars);
      %  res_mu_u(ub_vars) = mu.*ones(num_ub_vars,1)+dx(ub_vars).*dz_u(ub_vars)+(ub(ub_vars)-x(ub_vars)).*z_u(ub_vars);       
        res_mu_l(lb_vars) = mu.*ones(num_lb_vars,1)-dx(lb_vars).*dz_l(lb_vars);
        res_mu_u(ub_vars) = mu.*ones(num_ub_vars,1)+dx(ub_vars).*dz_u(ub_vars);       
        if (Newton_struct.la_mode == "exact")
            [dx_c,dy_c,dz_l_c,dz_u_c,instability] = ...
            Newton_struct.Newton_backsolve(NSdata,zeros(m,1),zeros(n,1),res_mu_l,res_mu_u);
        elseif (Newton_struct.la_mode == "inexact")
            [droptol] = Newton_struct.Precond_struct.set_out_precond_tol(droptol,iterlin,maxit_Krylov,nnzL,roof);
            PS.droptol = droptol; % Keep for monitoring purposes.
            [dx_c,dy_c,dz_l_c,dz_u_c,instability,iterlin,drop_direction,Reinforce_Inner_IR] = ...
            Newton_struct.Newton_itersolve(predictor,NSdata,PS,zeros(m,1),zeros(n,1), ...
                                           res_mu_l,res_mu_u,maxit_Krylov,Krylov_tol);
            % Keep track of the number of iterations performed. If an iter. is repeated, keep track of it.
            itertot(2*ctr_struct.IP_iter) = itertot(2*ctr_struct.IP_iter) + iterlin;  
            if (drop_direction == true)
                [ctr_struct,mu,delta,rho] = control_direction_accuracy(ctr_struct,mu,mu_prev,delta,rho,predictor);
                if (ctr_struct.opt == 5)    % 5 signifies failure due to inaccuracy
                    break;
                end
                continue;
            end
        end      
        if (instability == true)        % Checking if the matrix is too ill-conditioned. Mitigate it.
            [ctr_struct,delta,rho,reg_limit] = control_ill_conditioning(ctr_struct,rho,delta,reg_limit,tol,predictor);
            if (ctr_struct.opt == 4)    % 4 signifies failure due to numerical inaccuracy.
                break;
            end
            continue;
        end
        ctr_struct.retry_c = 0;
        if (Newton_struct.la_mode == "inexact" && Reinforce_Inner_IR)
            ctr_struct.set_proximal_inner_IR = true;
        end
        dx = dx + dx_c;
        dy = dy + dy_c;
        dz_l = dz_l + dz_l_c;
        dz_u = dz_u + dz_u_c;
        % ____________________________________________________________________________________________________________ %

        % ============================================================================================================ %
        % Step in the lower-or-upper (primal) and non-negativity orthant (dual).
        % ------------------------------------------------------------------------------------------------------------ %
        [alpha_x,alpha_z] = step_to_the_boundary(x,z_l,z_u,dx,dz_l,dz_u,lb,ub,lb_vars,ub_vars);
        % ____________________________________________________________________________________________________________ %
    end
    % ================================================================================================================ %
    % Make the step and compute the approximate rate of decrease (or increase) or mu.
    % ---------------------------------------------------------------------------------------------------------------- %
    x = x+alpha_x.*dx; y = y+alpha_z.*dy; z_l = z_l+alpha_z.*dz_l; z_u = z_u+alpha_z.*dz_u;
    if (num_lb_vars + num_ub_vars > 0)  % Only if we have inequality constraints.
        mu_prev = mu;
        mu = (x(lb_vars)-lb(lb_vars))'*z_l(lb_vars) + (ub(ub_vars)-x(ub_vars))'*z_u(ub_vars);
        %mu = max(mu/(num_lb_vars + num_ub_vars),0.1*mu_prev);
        mu = mu/(num_lb_vars + num_ub_vars);
        mu_rate = min(abs((mu-mu_prev)/max(mu,mu_prev)),0.9); 
        mu_rate = max(mu_rate,0.2);
    else
        mu_rate = 0.9;      % Arbitrary choice.
    end
    % ________________________________________________________________________________________________________________ %
    
    % ================================================================================================================ %
    % Computing the new non-regularized residuals. If the overall error is decreased, for the primal and dual 
    % residuals, we accept the new estimates for the Lagrange multipliers and primal optimal solution respectively.
    % If not, we keep the estimates constant. However, 
    % we continue decreasing the penalty parameters, limiting the decrease to the value of the minimum pivot
    % of the LDL^T decomposition (to ensure single pivots).
    % ---------------------------------------------------------------------------------------------------------------- %
    new_nr_res_p = b-A_struct.A(x);
    new_nr_res_d = c + Q_struct.Q(x) - A_struct.A_tr(y) - z_l + z_u;
    PMM_struct = build_PMM_structure(zeta,lambda,delta,rho,reg_limit,mu_rate);
    [PMM_struct,ctr_struct] = update_PMM_parameters(ctr_struct,nr_res_p,new_nr_res_p,nr_res_d,new_nr_res_d,x,y,PMM_struct);
    lambda = PMM_struct.lambda;     zeta = PMM_struct.zeta;     
    delta = PMM_struct.delta;       rho = PMM_struct.rho; 
    % ________________________________________________________________________________________________________________ %

    % ================================================================================================================ %
    % Print iteration output.  
    % ---------------------------------------------------------------------------------------------------------------- %
    pres_inf = norm(new_nr_res_p);
    dres_inf = norm(new_nr_res_d);  
    if (Newton_struct.la_mode == "exact")
        IP_PMM_output(fid,printlevel,ctr_struct.IP_iter,pres_inf,dres_inf,mu,alpha_x,alpha_z,delta,rho);
    elseif (Newton_struct.la_mode == "inexact")
        Krylov_its = Krylov_its + itertot(2*ctr_struct.IP_iter-1) + itertot(2*ctr_struct.IP_iter);
        IP_PMM_output(fid,printlevel,ctr_struct.IP_iter,pres_inf,dres_inf,mu,alpha_x,alpha_z,delta,rho,...
                      Newton_struct.la_mode,Krylov_its,nnzL);
    end
    % ________________________________________________________________________________________________________________ %
end % while (IPiter < maxit)

% The IPM has terminated because the solution accuracy is reached or the maximum number 
% of iterations is exceeded, or the problem under consideration is infeasible. Print result. 
IP_iter = ctr_struct.IP_iter;
if (Newton_struct.la_mode == "exact") Krylov_its = IP_iter; end
opt =  ctr_struct.opt;
fprintf('iterations: %4d\n', IP_iter);
fprintf('primal feasibility: %8.2e\n', norm(A_struct.A(x)-b));
fprintf('dual feasibility: %8.2e\n', norm(A_struct.A_tr(y)+z_l-z_u-c-Q_struct.Q(x)));
fprintf('complementarity: %8.2e\n', mu);  
end
% ******************************************************************************************************************** %
% END OF FILE
% ******************************************************************************************************************** %