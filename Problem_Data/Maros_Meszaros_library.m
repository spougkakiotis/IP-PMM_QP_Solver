function [solution_statistics] = Maros_Meszaros_library(set,type,tol,max_IP_iter,printlevel,la_mode,fid)
% ==================================================================================================================== %
% This function loads various QP problems from the Maros Meszaros library and solves them using IP_PMM.
%
% INPUT: 
% The first input argument---set---must contain a set of integer numbers 
% between 1 and 122 (number of Netlib problems). It must be given 
% as a row vector! It is then used 
% to call the respective problems from the test set. For example, for set = {1,2}
% the function attempts to solve only the first two problems of the test set.
% This argument cannot be empty.
%
% The second input argument---type---must be empty, or take between two values:
%               1) "standard",
%               2) "CONT",
% specifying whether we want to solve the standard or the CONT (5 problems) test set.
% (Default at "standard")
% The third input argument---tol---indicates the tolerance required.
% (Default at 1e-6)
% The fourth argument---max_IP_iter---specifies the maximum allowed number of IP iterations.
% (Default at 100)
% The fifth argument---printlevel---specifies the printing options. See the documentation 
% of IP-PMM for more. (Default at 1).
% The sixth argument---la_mode---specifies whether we want to use "exact" linear algebra
% (i.e. factorization) or "inexact" (i.e. Krylov subspace methods).
% Finally, the last argument---fid---specifies the file on which 
% the algorithmic printing is done.
%
% OUTPUT: The output is given in the form of a struct, collecting various statistics 
%         from the run of the method, such as:
%               total_time  -> Time needed to solve all specified problems
%               total_iters -> Number of IP-PMM iterations performed
%               problems_converged -> Number of problems converged
%               tol -> the tolerance used
%               problem_attempts -> Number of problems attempted to be solved
%               success_rate -> problems_converged/problem_attempts
%               total_Krylov_iters -> total number of Krylov iterations needed
%               objective_values -> a vector containing the optimal objective values of all problems solved
%               status_vec -> a vector containing the termination status of IP-PMM for each problem
%               max_nnzL -> a vector containing the maximum nnz used for preconditioning for each problem
%               solution_struct -> the solution of the last problem solved
%               problem_names -> contains the names of the problems solved
% ____________________________________________________________________________________________________________________ %
    if (isempty(set) || nargin < 1)
        error('Maros/Meszaros problem(s) not specified.\n');
    end
    if (nargin < 2 || isempty(type))        type = "standard"; end
    if (nargin < 3 || isempty(tol))         tol = 1e-6;        end
    if (nargin < 4 || isempty(max_IP_iter)) max_IP_iter = 100; end
    if (nargin < 5 || isempty(printlevel))  printlevel = 1;    end

    %The path on which all the netlib problems lie
    if (type == "standard")
        QP_problems_path = 'Problem_Data/Maros_Meszaros_QP_Collection/maros'; 
    elseif (type == "CONT")
        QP_problems_path = 'Problem_Data/Maros_Meszaros_QP_Collection/maros_CONT'; 
    else
        error("Incorrect input argument for type. Use either %s or %s.\n","standard","CONT");
    end

    %Finds all the QP problems and stores their names in a struct
    d = dir(fullfile(QP_problems_path,'*.mat')); 
    %Each indice i=1..num_of_QP_files gives the name of each QP problem though d(i).name
    
    %Open the file to write the results
    fileID = fopen('./Output_files/QP_problems_tabular_fortmat_final_results.txt','a+');
    
    solution_statistics = struct();
    solution_statistics.total_iters = 0; solution_statistics.total_time = 0;
    solution_statistics.problems_converged = 0; solution_statistics.tol = tol;
    solution_statistics.problem_attempts = 0; solution_statistics.success_rate = 0;
    solution_statistics.total_Krylov_iters = 0;
    solution_statistics.objective_values = zeros(size(set));
    solution_statistics.status_vec = zeros(size(set));
    solution_statistics.problem_names = strings(size(set));     % To keep the names of the problems solved.
    solution_statistics.max_nnzL = zeros(size(set));            % To keep the maximum nnz of a Cholesky factor found.
    scaling_direction = 'l'; scaling_option = 3; 
    
    model = struct();
    if (type == "standard") % b is not present but will be added
        fields = ["H","A","g","xl","xu","al","au","g0","name","b"];
    else % g0, name and b are not present but will be added
        fields = ["Q","A","c","lb","ub","rl","ru","g0","name","b"];
    end
    it_counter = 0;
    for k = set
        it_counter = it_counter + 1;
        if (isfield(model,fields)) %If any of the fields is missing, dont remove anything
            model = rmfield(model,fields); %Remove all fields before loading new ones
        end
        model = load(fullfile(QP_problems_path,d(k).name));
        if (type == "CONT")
            model.name = d(k).name;
            model.g0 = 0;
        end
    
        % Convert to problem to a box form which coincides with IP-PMM's input.
        [model.(fields(2)),model.(fields(1)),model.(fields(10)),...
        model.(fields(3)),model.(fields(4)),model.(fields(5))] =  ...
        Maros_Meszaros_Convert_to_Box_Form(model.(fields(2)),model.(fields(1)),model.(fields(3)), ...
                         model.(fields(6)), model.(fields(7)), model.(fields(4)), model.(fields(5)));
                     
        D = Scale_the_problem(model.(fields(2)),scaling_option,scaling_direction);
        time = 0;  tic;
        [solution_struct] = Set_Up_IP_PMM(model.(fields(2)),model.(fields(1)),model.(fields(10)),model.(fields(3)),...
                                       model.(fields(8)),model.(fields(4)),model.(fields(5)),tol,max_IP_iter,printlevel,la_mode,fid);                                      
        time = time + toc;
        solution_statistics.total_iters = solution_statistics.total_iters + solution_struct.IP_iter;
        solution_statistics.total_time = solution_statistics.total_time + time;  
        solution_statistics.total_Krylov_iters = solution_statistics.total_Krylov_iters + solution_struct.Krylov_its;
        solution_statistics.max_nnzL(it_counter) = solution_struct.max_nnzL;
        if (solution_struct.opt == 1)                                       % Success
           solution_statistics.problems_converged = solution_statistics.problems_converged + 1;
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
        fprintf(fileID,'%s & %d & %d & opt = %s  \n',model.(fields(9)), solution_struct.IP_iter,...
                                                     time, string(solution_struct.opt == 1)); 
        solution_statistics.objective_values(it_counter) = solution_struct.obj_val;
        solution_statistics.status_vec(it_counter) = solution_struct.opt;
        solution_statistics.problem_names(it_counter) = string(model.(fields(9)));
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

