function [c, A, b, lb, ub, obj_const_term] = LP_Convert_to_IP_PMM_Form(c, A, b, lb, ub, sense)
% ==================================================================================================================== %
% LP_Convert_to_IP_PMM_Form( c, A, b, lb, ub, sense ):
% This function takes as input the data of an LP in the following form:
%                       min    c^T x
%                       s.t.   Ax {<=,=,>=} rhs,
%                              lb <= x <= ub
% and transforms it in a box form, that is:
%                       min    c_bar^T x_bar
%                       s.t.   A_bar x_bar = b
%                              lb_bar <= x_bar <= ub_bar,
% where x_bar has been augmented by certain slack variables, and any variable has either a lower or an upper bound.
% The objective value does not change by this reformulation.
%
% Author: Spyridon Pougkakiotis.
% ____________________________________________________________________________________________________________________ %
    obj_const_term = 0;         % Any constant term of the objective removed due to the reformulation.
    % ================================================================================================================ %
    % Test input data, dimensions, e.t.c.
    % ---------------------------------------------------------------------------------------------------------------- %
    [m,n] = size(A);
    if (size(lb,2) > 1)
        lb = lb';
    elseif (size(ub,2) > 1)
        ub = ub';
    elseif (size(sense,2) > 1)
        sense = sense';
    elseif (size(c,2) > 1)
        c = c';
    elseif (~issparse(A))
        A = sparse(A);
    end

    if (size(c,1) ~= n || size(sense,1) ~= m || size(lb,1) ~= n || size(ub,1) ~= n)
        error("Incorrect input dimensions")
    end
    % ________________________________________________________________________________________________________________ %


    % ================================================================================================================ %
    % Initialization.
    % ---------------------------------------------------------------------------------------------------------------- %
    num_of_slacks = 0; % Counter for the slack variables to be added in the inequality constraints.
    [rows,cols,v] = find(A);
    % ________________________________________________________________________________________________________________ %

    % ================================================================================================================ %
    % Make all the constraints to be of equality type (add slack variables)
    % ---------------------------------------------------------------------------------------------------------------- %
    for i = 1:m    
        if ( sense(i) == '<') % add a slack of the form +x_slack.               
           num_of_slacks = num_of_slacks + 1;
           rows = [rows; i];
           cols = [cols; n + num_of_slacks];
           v = [v; 1]; % assign 1 in the element A(i,n+num_of_slacks).       
           ub = [ub; Inf];
           lb = [lb; 0];
        elseif ( sense(i) == '>') % add a slack of the form -x_slack.                 
           num_of_slacks = num_of_slacks + 1;
           rows = [rows; i];
           cols = [cols;  n + num_of_slacks];
           v = [v; -1]; %assign -1 in the element A(i,n+num_of_slacks).    
           ub = [ub; Inf];
           lb = [lb; 0];
        end   
    end
    % ________________________________________________________________________________________________________________ %
    c = [c; zeros(num_of_slacks,1)];
    A = sparse(rows,cols,v,m,n+num_of_slacks);
    % ================================================================================================================ %
    % Find all variables that are constant and remove them from the problem
    % ---------------------------------------------------------------------------------------------------------------- %
    constant_vars = (lb == ub);  % Index of variables that are constant.
    num_constant_vars = nnz(constant_vars);
    if (num_constant_vars > 0)
        non_constant_vars = ~constant_vars; 
        obj_const_term = c(constant_vars)'*lb(constant_vars);
        b = b - A(:,constant_vars)*lb(constant_vars);
        c = c(non_constant_vars);
        lb = lb(non_constant_vars);
        ub = ub(non_constant_vars);
        A = A(:,non_constant_vars);
    end
    % ________________________________________________________________________________________________________________ %
    % ================================================================================================================ %
    % Introduce constraints so that no variable has both upper and lower bounds.
    % ---------------------------------------------------------------------------------------------------------------- %
    [rows,cols,v] = find(A);
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
    % ________________________________________________________________________________________________________________ %
    
    
end

