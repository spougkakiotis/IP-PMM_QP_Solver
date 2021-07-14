function [A, Q, b, c, lb, ub] = Maros_Meszaros_Convert_to_Box_Form(A, Q, c, al, au, lb, ub)
% ==================================================================================================================== %
% Maros_Meszaros_Convert_to_Box_Form(A, Q, c, al, au, lb, ub):
% This function takes as input the data of a QP problem from the Maros Meszaros library in the following form:
%                       min    c^T x + (1/2)x^T Q x
%                       s.t.   al <= Ax <= au,
%                              lb <= x <= ub
% and transforms it in a semi-standard form, that is:
%                       min    c_bar^T x + (1/2)x^T Q_bar x
%                       s.t.   A_bar x = b
%                              lb_bar <= x_bar <= ub_bar
%
% This does not alter the objective value of the problem.
%
% Author: Spyridon Pougkakiotis
% ==================================================================================================================== %

    % ================================================================================================================ %
    % Test input data, dimensions, e.t.c.
    % ---------------------------------------------------------------------------------------------------------------- %
    [m, n] = size(A);
    if (size(lb,2) > 1) lb = lb';       end
    if (size(ub,2) > 1) ub = ub';       end
    if (size(au,2) > 1) au = au';       end
    if (size(al,2) > 1) al = al';       end
    if (size(c,2) > 1)  c = c';         end
    if (~issparse(A))   A = sparse(A);  end
    if (~issparse(Q))   Q = sparse(Q);  end


    if (size(c,1) ~= n || size(al,1) ~= m || size(au,1) ~= m || size(lb,1) ~= n         ...
                       || size(ub,1) ~= n || size(Q,1) ~= n || size(Q,2) ~= n)
        error("Incorrect input dimensions.\n")
    end
    % ________________________________________________________________________________________________________________ %

    % ================================================================================================================ %
    % Initialization.
    % ---------------------------------------------------------------------------------------------------------------- %
    num_of_slacks = 0;        % Counter for the slack variables to be added in the inequality constraints.
    extra_constraints = 0;    % Counter for the extra constraints to be added in case of double bounds.
    [rows,cols,v] = find(A);
    if (size(rows,2) > 1)
        rows = rows';
    end
    if (size(cols,2) > 1)
        cols = cols';
    end
    if (size(v,2) > 1)
        v = v';
    end
    b = zeros(m,1);
    % ________________________________________________________________________________________________________________ %

    % ================================================================================================================ %
    % Make all the constraints to be of equality type (add slack variables)
    % ---------------------------------------------------------------------------------------------------------------- %
    for i = 1:m   
        if (au(i) ~= al(i))
            if (au(i) == Inf)
                % we only have a lower bound, and hence we should add a slack of the form -x_slack
                num_of_slacks = num_of_slacks + 1;
                rows = [rows; i];
                cols = [cols; n + num_of_slacks];
                v = [v; -1];        % assign -1 in the element A(i,n+num_of_slacks) 
                b(i) = al(i);       % Fix the RHS    
                lb = [lb; 0];
                ub = [ub; Inf];
            elseif (al(i) == -Inf)
                % we only have an upper bound, and hence we should add a slack of the form +x_slack
                num_of_slacks = num_of_slacks + 1;
                rows = [rows; i];
                cols = [cols; n + num_of_slacks];
                v = [v; 1];         % assign 1 in the element A(i,n+num_of_slacks)
                b(i) = au(i);       % Fix the RHS     
                lb = [lb; 0];
                ub = [ub; Inf];
            else
                % transform al <=Axi <=au to Axi' = aui, Axi' = ali
                extra_constraints = extra_constraints + 1;
                k_max = size(cols,1);
                for k = 1:k_max
                    if (rows(k) == i)
                        cols = [cols; cols(k)];
                        rows = [rows; m + extra_constraints];
                        v = [v; v(k)];
                    end
                end

                % treat the case of the upper bound
                num_of_slacks = num_of_slacks + 1;
                rows = [rows; i];
                cols = [cols; n + num_of_slacks];
                v = [v; 1];         % assign 1 in the element A(i,n+num_of_slacks)
                b(i) = au(i); % Fix the RHS

                % Now add a new constraint that will treat the case of the LB
                num_of_slacks = num_of_slacks + 1;
                rows = [rows; m + extra_constraints];
                cols = [cols; n + num_of_slacks];
                v = [v; -1];
                b = [b; al(i)]; % The RHS of the extra constraint     
                lb = [lb; 0; 0];
                ub = [ub; Inf; Inf];
            end      
        else
            b(i) = al(i); % Already an equality constraint.
        end

    end
    % ________________________________________________________________________________________________________________ %
    c = [c; zeros(num_of_slacks,1)];
    A = sparse(rows,cols,v,m + extra_constraints, n + num_of_slacks); % Renew the matrix to incude new constraints.
    Q = [Q sparse(n,num_of_slacks); sparse(num_of_slacks,n) sparse(num_of_slacks,num_of_slacks)];
end

