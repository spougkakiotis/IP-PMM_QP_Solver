function [solution_statistics] = Netlib_library(set,type,tol,max_IP_iter,printlevel)
% ==================================================================================================================== %
% This function loads various NETLIB problems and solves them using IP_PMM.
%
% INPUT: 
% The first input argument---set---must contain a set of integer numbers 
% between 1 and 96 (number of Netlib problems). It must be given 
% as a row vector! It is then used 
% to call the respective problems from the test set. For example, for set = {1,2}
% the function attempts to solve only the first two problems of the test set.
% This argument cannot be empty.
%
% The second input argument---type---must be empty, or take between two values:
%               1) "standard",
%               2) "presolved",
% specifying whether we want to solve the standard or the presolved test set.
% (Default at "standard")
% The third input argument---tol---indicates the tolerance required.
% (Default at 1e-6)
% The fourth argument---max_IP_iter---specifies the maximum allowed number of IP iterations.
% (Default at 100)
% The fifth argument---printlevel---specifies the printing options. See the documentation 
% of IP-PMM for more. (Default at 1).
%
% OUTPUT: The output is given in the form of a struct, collecting various statistics 
%         from the run of the method, such as:
%               total_time  -> Time needed to solve all specified problems
%               total_iters -> Number of IP-PMM iterations performed
%               problems_converged -> Number of problems converged
%               problem_attempts -> Number of problems attempted to be solved
%               success_rate -> problems_converged/problem_attempts
%               solution_struct -> the solution of the last problem solved
% ____________________________________________________________________________________________________________________ %
    if (isempty(set) || nargin < 1)
        error('Netlib problem(s) not specified.\n');
    end
    if (nargin < 2 || isempty(type))        type = "standard"; end
    if (nargin < 3 || isempty(tol))         tol = 1e-6;        end
    if (nargin < 4 || isempty(max_IP_iter)) max_IP_iter = 100; end
    if (nargin < 5 || isempty(printlevel))  printlevel = 1;    end
    %The path on which all the netlib problems lie
    Netlib_path = './Problem_Data/Netlib_LP_Collection'; 
    %Finds all the Netlib problems and stores their names in a struct
    d = dir(fullfile(Netlib_path,'*.mat'));
    %Open the file to write the results
    fileID = fopen('./Output_files/Netlib_tabular_format_IP_PMM_runs.txt','a+');

    solution_statistics = struct();
    solution_statistics.total_iters = 0; solution_statistics.total_time = 0;
    solution_statistics.problems_converged = 0; solution_statistics.tol = tol;
    solution_statistics.problem_attempts = 0; solution_statistics.success_rate = 0;
    
    scaling_direction = 'l'; scaling_mode = 3; 
    for k = set %Each indice k=1..num_of_netlib_files gives the name of each netlib problem through d(i).name
        load(fullfile(Netlib_path,d(k).name));                              % Loads two structs: The standard model and 
                                                                            % a presolved model.
        %fields = {'A','obj','sense','rhs','lb','ub','vtype','modelname','varnames','constrnames'};
        if (type == "standard")
            problem = model;
        elseif (type == "presolved")
            problem = presolved_model;
        end
        [problem.obj, problem.A, problem.rhs, problem.lb, problem.ub, obj_const_term] = ...
        LP_Convert_to_IP_PMM_Form(problem.obj, problem.A,  problem.rhs, problem.lb, problem.ub, problem.sense);
        n = size(problem.lb,1);
        m = size(problem.rhs,1);
        Q = sparse(n,n);
        if (scaling_direction == 'l')
            [D,~] = Scale_the_problem(problem.A,scaling_mode,scaling_direction);
            problem.A = spdiags(D,0,m,m)*problem.A;                         % Apply the left scaling.
            problem.rhs = problem.rhs.*D;
        end
        time = 0;  tic;
        [solution_struct] = Set_Up_IP_PMM(problem.A,Q,problem.rhs,problem.obj,obj_const_term,...
                                          problem.lb,problem.ub,tol,max_IP_iter,printlevel);
        time = time + toc;
        solution_statistics.total_iters = solution_statistics.total_iters + solution_struct.IP_iter;
        solution_statistics.total_time = solution_statistics.total_time + time;  
        
        if (solution_struct.opt == 1)                                       % Success
           solution_statistics.problems_converged = solution_statistics.problems_converged + 1;
           fprintf(fileID,'The optimal solution objective is %d.\n',solution_struct.obj_val);
        elseif (solution_struct.opt == 2 || solution_struct.opt == 3)       % Infeasibility
           fprintf(fileID,'Infeasibility has been detected.\n');
        else                                                                % Reached maximum iterations
           fprintf(fileID,'Maximum number of iterations reached.\n');
        end
        fprintf(fileID,'%s & %d & %d & opt = %s  \n',problem.modelname, solution_struct.IP_iter,...
                                                     time, string(solution_struct.opt == 1)); 
        
        solution_statistics.problem_attempts = solution_statistics.problem_attempts + 1;
        solution_statistics.success_rate = ...
            solution_statistics.problems_converged/solution_statistics.problem_attempts;
    end
    solution_statistics.solution_struct= solution_struct;
    fprintf(fileID,['The total iterates were: %d and the total time was %d.',...
                    '%d problems converged.\n'],solution_statistics.total_iters,...
                    solution_statistics.total_time,solution_statistics.problems_converged);
    fclose(fileID);
end

