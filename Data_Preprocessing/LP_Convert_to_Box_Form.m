function [c, A, b, lb, ub] = LP_Convert_to_Box_Form(c, A, b, lb, ub, sense)
% ==================================================================================================================== %
% LP_Convert_to_Box_Form(c, A, b, lb, ub, sense):
% This function takes as input the data of a Netlib LP in the following form:
%                       min    c^T x
%                       s.t.   Ax {<=,=,>=} rhs,
%                              lb <= x <= ub
% and transforms it in a box form, that is:
%                       min    c_bar^T x_bar
%                       s.t.   A_bar x_bar = b
%                              lb_bar <= x_bar <= ub_bar,
% where x_bar has been augmented by certain slack variables.
% The objective value does not change by this reformulation.
%
% Author: Spyridon Pougkakiotis.
% ____________________________________________________________________________________________________________________ %
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
    if (size(rows,1) == 1)
        rows = rows';
    end
    if (size(cols,1) == 1)
        cols = cols';
    end
    if (size(v,1) == 1)
        v = v';
    end
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
end

