function [A, Q, b, c, lb, ub, obj_const_term] = QP_Convert_to_L_Or_U_Form(A, Q, b, c, lb, ub)
% ==================================================================================================================== %
% QP_Convert_to_L_Or_U_Form( c, A, b, lb, ub, sense ):
% This function takes as input the data of an LP in the following form:
%                       min    c^T x + (1/2) x^T Q x
%                       s.t.   Ax = b,
%                              lb <= x <= ub
% and transforms it in a box form, that is:
%                       min    c_bar^T x_bar
%                       s.t.   A_bar x_bar = b
%                              x_bar {>=, <=} {lb_bar,ub_bar}
% where x_bar has been augmented by certain slack variables, and any variable i has either a lower (lb_bar(i))
% or an upper bound (ub_bar(i)). In particular, 
%      if (lb_bar(i) >= -Inf) then ub_bar(i) = Inf
%  and if (ub_bar(i) < Inf)   then lb_bar(i) = -Inf.
% Any constant variable is deleted and the method returns a "obj_const_term" which is the constant 
% term that needs to be added to the objective.
%
% Author: Spyridon Pougkakiotis.
% ____________________________________________________________________________________________________________________ %
    obj_const_term = 0;         % Any constant term of the objective removed due to the reformulation.
    % ================================================================================================================ %
    % Test input data, dimensions, e.t.c.
    % ---------------------------------------------------------------------------------------------------------------- %
    [m,n] = size(A);
    if (size(lb,2) > 1)     lb = lb';       end
    if (size(ub,2) > 1)     ub = ub';       end
    if (size(c,2) > 1)      c = c';         end
    if (~issparse(A))       A = sparse(A);  end
    if (~issparse(Q))       Q = sparse(Q);  end

    if (size(c,1) ~= n || size(b,1) ~= m || size(Q,1) ~= n || size(Q,2) ~= n || size(lb,1) ~= n || size(ub,1) ~= n)
        error("Incorrect input dimensions")
    end
    % ________________________________________________________________________________________________________________ %

    % ================================================================================================================ %
    % Find all variables that are constant and remove them from the problem
    % ---------------------------------------------------------------------------------------------------------------- %
    constant_vars = (lb == ub);  % Index of variables that are constant.
    num_constant_vars = nnz(constant_vars);
    if (num_constant_vars > 0)
        non_constant_vars = ~constant_vars; 
        obj_const_term = c(constant_vars)'*lb(constant_vars) + ...
                       (1/2)*((lb(constant_vars))'*(Q(constant_vars,constant_vars)*lb(constant_vars)));
        b = b - A(:,constant_vars)*lb(constant_vars);
        c = c(non_constant_vars) + Q(non_constant_vars,constant_vars)*lb(constant_vars);
        lb = lb(non_constant_vars);
        ub = ub(non_constant_vars);
        A = A(:,non_constant_vars);
        Q = Q(non_constant_vars,non_constant_vars);
    end
    % ________________________________________________________________________________________________________________ %

    % ================================================================================================================ %
    % Introduce constraints so that no variable has both upper and lower bounds.
    % ---------------------------------------------------------------------------------------------------------------- %
    [rows,cols,v] = find(A);
    if (size(rows,1) == 1)
        rows = rows';
    end
    if (size(cols,1) == 1)
        cols = cols';
    end
    if (size(v,1) == 1)
        v = v';
    end
    num_extra_constraints = 0;
    n = size(A,2);
    for i = 1:n
        if ((lb(i) > -Inf) && (ub(i) < Inf))  % We have both upper and lower bound -> New constraint to treat UB
            num_extra_constraints = num_extra_constraints + 1;
            b = [b; ub(i)];  %The RHS of extra constraint x_i + w_i = ub_i - lb_i
            rows = [rows; m + num_extra_constraints ; m + num_extra_constraints];
            cols = [cols; i ; n + num_extra_constraints];
            v = [v; 1; 1];   % Assigns ones in the element A(m+extra_constr,i) and A(m+extra_constr,n+extra_constr) 
            ub(i) = Inf;     % Forget about the upper bound of the original variable
            lb = [lb; 0];    % Lower bound for the slack variable
            ub = [ub; Inf];  % Upper bound for the slack variable
        end 
    end
    c = [c; zeros(num_extra_constraints,1)];
    A = sparse(rows,cols,v,m+num_extra_constraints,n+num_extra_constraints);
    Q = [Q                               sparse(n,num_extra_constraints); 
         sparse(num_extra_constraints,n) sparse(num_extra_constraints,num_extra_constraints)];
    % ________________________________________________________________________________________________________________ %    
end

