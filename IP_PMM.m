function [x,y,z_l,z_u,opt,IP_iter] = IP_PMM(A_struct,Q_struct,Newton_struct,b,c,lb,ub,tol,maxit,printlevel,fid)
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
%                                              3: print intermediate Krylov iterates
% IP_PMM(A, Q, Newton_struct, b, c, lb, ub, tol, max_it, printlevel,fid): prints on a file with indicator fid.
% OUTPUT: [x,y,z,opt,iter], where:
%         x: primal solution
%         y: Lagrange multiplier vector
%         z_l: dual slack variables corresponding to the lower bound constraints
%         z_u: dual slack variables corresponding to the upper bound constraints
%         opt: true if problem was solved to optimality, false if problem not solved or found infeasible.
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
    A_struct.NE_diag = zero_function;
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
if (nargin < 6  || isempty(lb))            lb = -Inf.*ones(n,1); end
if (nargin < 7  || isempty(ub))            ub = Inf.*ones(n,1);  end
if (nargin < 8 || isempty(tol))            tol = 1e-4;           end
if (nargin < 9 || isempty(maxit))          maxit = 100;          end
if (nargin < 10 || isempty(printlevel))    printlevel = 1;       end
if (nargin < 11 || isempty(fid))           fid = 1;       end
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
%x = 1.*ones(n,1);
%y = 1.*ones(m,1);
%z_l = 1.*ones(n,1);
%z_u = zeros(n,1);
lambda = y;   zeta = x;                                   % Initial estimate of the primal and dual optimal solution.
if (num_lb_vars + num_ub_vars > 0)                        % Initial value of mu (when bound constraints are present).
    mu = (x(lb_vars)-lb(lb_vars))'*z_l(lb_vars) + ...
         (ub(ub_vars) - x(ub_vars))'*z_u(ub_vars);
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
IP_PMM_header(fid,printlevel);                            % Set the header for printing
delta = 8;  rho = 8;                                      % Initial dual and primal regularization value.
reg_limit = max(5*tol*(1/max(A_struct.A_1_norm*A_struct.A_inf_norm,...
                Q_struct.Q_1_norm*Q_struct.Q_inf_norm)),5*1e-10);
reg_limit = min(reg_limit,1e-8);                          % Controlled perturbation.
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
    % ================================================================================================================ %
    % Check terminatio criteria
    % ---------------------------------------------------------------------------------------------------------------- %
    if (norm(nr_res_p)/(1+norm(b)) < tol && norm(nr_res_d)/(1+norm(c)) < tol &&  mu < tol )
        fprintf('optimal solution found\n');
        ctr_struct.opt = 1;
        break;
    end
    if ((norm(y-lambda)> 10^10 && norm(res_p) < tol && no_dual_update > 5)) 
        fprintf('The primal-dual problem is infeasible\n');
        ctr_struct.opt = 2;
        break;
    end
    if ((norm(x-zeta)> 10^10 && norm(res_d) < tol && no_primal_update > 5))
        fprintf('The primal-dual problem is infeasible\n');
        ctr_struct.opt = 3;
        break;
    end
    % ________________________________________________________________________________________________________________ %
    ctr_struct.IP_iter = ctr_struct.IP_iter+1;
    % ================================================================================================================ %
    % Avoid the possibility of converging to a local minimum -> Decrease the minimum regularization value.
    % ---------------------------------------------------------------------------------------------------------------- %
  %  if (ctr_struct.no_primal_update > 5 && rho == reg_limit && reg_limit ~= 5*1e-13)
  %      reg_limit = 5*1e-13;
  %      ctr_struct.no_primal_update = 0;
  %      ctr_struct.no_dual_update = 0;
  %  elseif (ctr_struct.no_dual_update > 5 && delta == reg_limit && reg_limit ~= 5*1e-13)
  %      reg_limit = 5*1e-13;
  %      ctr_struct.no_primal_update = 0;
  %      ctr_struct.no_dual_update = 0;
  %  end
    % ________________________________________________________________________________________________________________ %
    % ================================================================================================================ %
    % Compute the Newton factorization.
    % ---------------------------------------------------------------------------------------------------------------- %
    pivot_thr = reg_limit/5;
    NS = Newton_struct.Newton_fact(x,z_l,z_u,delta,rho,lb,ub,lb_vars,ub_vars,pivot_thr);
    % ________________________________________________________________________________________________________________ %
 
    % ================================================================================================================ %
    % Predictor step: Solve the Newton system and compute a centrality measure.
    % ---------------------------------------------------------------------------------------------------------------- %
    if (num_lb_vars + num_ub_vars > 0)          % Predictor only in the case we have inequality constraints.
        res_mu_l(lb_vars) = (lb(lb_vars)-x(lb_vars)).*z_l(lb_vars);
        res_mu_u(ub_vars) = (x(ub_vars)-ub(ub_vars)).*z_u(ub_vars);
        % ============================================================================================================ %
        % Solve the Newton system with the predictor right hand side -> Optimistic view, solve as if you wanted to 
        %                                                               solve the original problem in 1 iteration.
        % ------------------------------------------------------------------------------------------------------------ %
        [dx,dy,dz_l,dz_u,instability] = Newton_backsolve(NS,res_p,res_d,res_mu_l,res_mu_u);
        if (instability == true) % Checking if the matrix is too ill-conditioned. Mitigate it.
            if (ctr_struct.retry_p < ctr_struct.max_tries)
                fprintf('The system is re-solved, due to bad conditioning  of predictor system.\n')
                delta = delta*100;  rho = rho*100;
                ctr_struct.IP_iter = ctr_struct.IP_iter -1;
                ctr_struct.retry_p = ctr_struct.retry_p + 1;
                reg_limit = min(reg_limit*10,tol);
                continue;
            else
                fprintf('The system matrix is too ill-conditioned.\n');
                ctr_struct.opt = 0;
                break;
            end
        end
        ctr_struct.retry_p = 0;
        % ____________________________________________________________________________________________________________ %
        
        % ============================================================================================================ %
        % Step in the box (primal) and non-negativity orthant (dual).
        % ------------------------------------------------------------------------------------------------------------ %
        pos_idx = false(n,1);   neg_idx = false(n,1);
        neg_idz_l = false(n,1); neg_idz_u = false(n,1);
        pos_idx(ub_vars) = dx(ub_vars) > 0;
        neg_idx(lb_vars) = dx(lb_vars) < 0;
        neg_idz_l(lb_vars) = dz_l(lb_vars) < 0;
        neg_idz_u(ub_vars) = dz_u(ub_vars) < 0;
        alphamax_x = min([1;(lb(neg_idx)-x(neg_idx))./dx(neg_idx);...
                         (ub(pos_idx)-x(pos_idx))./dx(pos_idx)]);
        alphamax_z = min([1;-z_l(neg_idz_l)./dz_l(neg_idz_l);...
                         -z_u(neg_idz_u)./dz_u(neg_idz_u)]);
        tau = 0.995;
        alpha_x = tau*alphamax_x;
        alpha_z = tau*alphamax_z;
        % ____________________________________________________________________________________________________________ %
        centrality_measure = (x(lb_vars)-lb(lb_vars) + alpha_x.*dx(lb_vars))'*(z_l(lb_vars) + alpha_z.*dz_l(lb_vars)) + ...
                             (ub(ub_vars) - x(ub_vars)-alpha_x.*dx(ub_vars))'*(z_u(ub_vars) + alpha_z.*dz_u(ub_vars));
        mu = (centrality_measure/((num_lb_vars + num_ub_vars)*mu))^3*mu; %(centrality_measure/(num_lb_vars + num_ub_vars));
    else
        dx = zeros(n,1); dy = zeros(m,1); dz_l = zeros(n,1); dz_u = zeros(n,1);
    end
    % ________________________________________________________________________________________________________________ %
        
    % ================================================================================================================ %
    % Corrector step: Solve Newton system with the corrector right hand side. Solve as if you wanted to direct the 
    %                 method in the center of the central path.
    % ---------------------------------------------------------------------------------------------------------------- %
        res_mu_l(lb_vars) = mu.*ones(num_lb_vars,1)-(x(lb_vars)-lb(lb_vars)).*z_l(lb_vars)-dx(lb_vars).*dz_l(lb_vars);
        res_mu_u(ub_vars) = mu.*ones(num_ub_vars,1)+(x(ub_vars)-ub(ub_vars)).*z_u(ub_vars)+dx(ub_vars).*dz_u(ub_vars);
        % ============================================================================================================ %
        % Solve the Newton system with the predictor right hand side -> Optimistic view, solve as if you wanted to 
        %                                                               solve the original problem in 1 iteration.
        % ------------------------------------------------------------------------------------------------------------ %
        [dx_c,dy_c,dz_l_c,dz_u_c,instability] = Newton_backsolve(NS,res_p,res_d,res_mu_l,res_mu_u);
        if (instability == true) % Checking if the matrix is too ill-conditioned. Mitigate it.
            if (ctr_struct.retry_c < ctr_struct.max_tries)
                fprintf('The system is re-solved, due to bad conditioning of corrector.\n')
                delta = delta*100;
                rho = rho*100;
                ctr_struct.IP_iter = ctr_struct.IP_iter -1;
                ctr_struct.retry_c = ctr_struct.retry_c + 1;
                mu = mu_prev;
                reg_limit = min(reg_limit*10,tol);
                continue;
            else
                fprintf('The system matrix is too ill-conditioned.\n');
                ctr_struct.opt = 0;
                break;
            end
        end
        ctr_struct.retry_c = 0;
        % ____________________________________________________________________________________________________________ %
        dx = dx + dx_c;
        dy = dy + dy_c;
        dz_l = dz_l + dz_l_c;
        dz_u = dz_u + dz_u_c;
    % ________________________________________________________________________________________________________________ %
    
    
    % ================================================================================================================ %
    % Compute the new iterate:
    % Determine primal and dual step length. Calculate "step to the boundary" alphamax_x and alphamax_z. 
    % Then choose 0 < tau < 1 heuristically, and set step length = tau * step to the boundary.
    % Step in the box (primal) and non-negativity orthant (dual).
    % ---------------------------------------------------------------------------------------------------------------- %
    if (num_lb_vars + num_ub_vars > 0)
        pos_idx = false(n,1);   neg_idx = false(n,1);
        neg_idz_l = false(n,1); neg_idz_u = false(n,1);
        pos_idx(ub_vars) = dx(ub_vars) > 0;
        neg_idx(lb_vars) = dx(lb_vars) < 0;
        neg_idz_l(lb_vars) = dz_l(lb_vars) < 0;
        neg_idz_u(ub_vars) = dz_u(ub_vars) < 0;
        alphamax_x = min([1;(lb(neg_idx)-x(neg_idx))./dx(neg_idx);...
                         (ub(pos_idx)-x(pos_idx))./dx(pos_idx)]);
        alphamax_z = min([1;-z_l(neg_idz_l)./dz_l(neg_idz_l);...
                         -z_u(neg_idz_u)./dz_u(neg_idz_u)]);        
        tau = 0.995;
        alpha_x = tau*alphamax_x;
        alpha_z = tau*alphamax_z;
    else
        alpha_x = 1;         % If we have no inequality constraints, Newton method is exact -> Take full step.
        alpha_z = 1;        
    end
    % ________________________________________________________________________________________________________________ %
    
    % ================================================================================================================ %
    % Make the step.
    % ---------------------------------------------------------------------------------------------------------------- %
    x = x+alpha_x.*dx; y = y+alpha_z.*dy; z_l = z_l+alpha_z.*dz_l; z_u = z_u+alpha_z.*dz_u;
    if (num_lb_vars + num_ub_vars > 0)  % Only if we have inequality constraints.
        mu_prev = mu;
        mu = (x(lb_vars)-lb(lb_vars))'*z_l(lb_vars) + (ub(ub_vars)-x(ub_vars))'*z_u(ub_vars);
        mu = mu/(num_lb_vars + num_ub_vars);
        mu_rate = max(abs((mu-mu_prev)/max(mu,mu_prev)),0.95);
        mu_rate = max(mu_rate,0.1);
    else
        mu_rate = 0.9;
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
    cond = 0.95*norm(nr_res_p) > norm(new_nr_res_p);
    if (cond)
        lambda = y;
        delta = max(reg_limit,delta*(1-mu_rate));  
    else
        delta = max(reg_limit,delta*(1-0.666*mu_rate));     % Slower rate of decrease, to avoid losing centrality.       
        ctr_struct.no_dual_update = ctr_struct.no_dual_update + 1;
    end
    cond = 0.95*norm(nr_res_d) > norm(new_nr_res_d);
    if (cond)
        zeta = x;
        rho = max(reg_limit,rho*(1-mu_rate));  
    else
        rho = max(reg_limit,rho*(1-0.666*mu_rate));         % Slower rate of decrease, to avoid losing centrality.     
        ctr_struct.no_primal_update = ctr_struct.no_primal_update + 1;
    end
    % ________________________________________________________________________________________________________________ %

    % ================================================================================================================ %
    % Print iteration output.  
    % ---------------------------------------------------------------------------------------------------------------- %
    pres_inf = norm(new_nr_res_p);
    dres_inf = norm(new_nr_res_d);  
    IP_PMM_output(fid,printlevel,ctr_struct.IP_iter,pres_inf,dres_inf,mu,alpha_x,alpha_z);
    % ________________________________________________________________________________________________________________ %
end % while (iter < maxit)

% The IPM has terminated because the solution accuracy is reached or the maximum number 
% of iterations is exceeded, or the problem under consideration is infeasible. Print result. 
IP_iter = ctr_struct.IP_iter;
opt =  ctr_struct.opt;
fprintf('iterations: %4d\n', IP_iter);
fprintf('primal feasibility: %8.2e\n', norm(A_struct.A(x)-b));
fprintf('dual feasibility: %8.2e\n', norm(A_struct.A_tr(y)+z_l-z_u-c-Q_struct.Q(x)));
fprintf('complementarity: %8.2e\n', mu);  
end
% ******************************************************************************************************************** %
% END OF FILE
% ******************************************************************************************************************** %