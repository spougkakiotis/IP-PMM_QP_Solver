function [solution_statistics] = SuiteSparse_library(set,tol,max_IP_iter,printlevel,la_mode,fid)
% ==================================================================================================================== %
% This function loads various NETLIB problems and solves them using IP_PMM.
%
% INPUT: 
% The first input argument---set---must contain a set of integer numbers 
% between 1 and 29 (number of Netlib problems). It must be given 
% as a row vector! It is then used 
% to call the respective problems from the test set. For example, for set = {1,2}
% the function attempts to solve only the first two problems of the test set.
% This argument cannot be empty.
%
% The second input argument---tol---indicates the tolerance required.
% (Default at 1e-6)
% The third argument---max_IP_iter---specifies the maximum allowed number of IP iterations.
% (Default at 100)
% The fourth argument---printlevel---specifies the printing options. See the documentation 
% of IP-PMM for more. (Default at 1).
% The fifth argument---la_mode---specifies whether we want to use "exact" linear algebra
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
        error('SuiteSparse problem(s) not specified.\n');
    end
    if (nnz(set) > 29)
        error('Too many problems requested. The problem set is not so large.\n');
    end
    if (nargin < 2 || isempty(tol))         tol = 1e-6;        end
    if (nargin < 3 || isempty(max_IP_iter)) max_IP_iter = 100; end
    if (nargin < 4 || isempty(printlevel))  printlevel = 1;    end
    if (nargin < 5 || isempty(la_mode))     la_mode = "exact"; end
    if (nargin < 6 || isempty(fid))         fid = 1;           end

    %The path on which all the netlib problems lie
    SuiteSparse_path = './Problem_Data/SuiteSparse_LP_Collection'; 
    %Finds all the Netlib problems and stores their names in a struct
    d = dir(fullfile(SuiteSparse_path,'*.mat'));
    %Open the file to write the results
    fileID = fopen('./Output_files/SuiteSparse_tabular_format_IP_PMM_runs.txt','a+');
    
    solution_statistics = struct();
    solution_statistics.total_iters = 0; solution_statistics.total_time = 0;
    solution_statistics.problems_converged = 0; solution_statistics.tol = tol;
    solution_statistics.problem_attempts = 0; solution_statistics.success_rate = 0;
    solution_statistics.total_Krylov_iters = 0;
    solution_statistics.objective_values = zeros(size(set));    % To keep objective values.
    solution_statistics.status_vec = zeros(size(set));          % To keep convergence status.
    solution_statistics.max_nnzL = zeros(size(set));            % To keep the maximum nnz of a Cholesky factor found.
    solution_statistics.problem_names = strings(size(set));     % To keep the names of the problems solved.

    it_counter = 0;
    %for k = set
    %    load(fullfile(SuiteSparse_path,d(k).name));        % Loads the standard model
    %    disp(Problem.name);
    %end
    for k = set %Each indice k=1..num_of_suitesparse_files gives the name of each netlib problem through d(i).name
        it_counter = it_counter + 1;
        load(fullfile(SuiteSparse_path,d(k).name));        % Loads the standard model
        obj_const_term = Problem.aux.z0;                   % Constant term to add to the objective.
        Problem.sense = zeros(size(Problem.A,1),1);         
        Problem.sense(:) = '=';                            % We only have equality constraints.
        [Problem.aux.c, Problem.A, Problem.b, Problem.aux.lo, Problem.aux.hi] = ...
        LP_Convert_to_Box_Form(Problem.aux.c, Problem.A,  Problem.b, Problem.aux.lo, Problem.aux.hi, Problem.sense);
        n = size(Problem.aux.lo,1);
        m = size(Problem.b,1);
        Q = sparse(n,n);
        [solution_struct] = Set_Up_IP_PMM(Problem.A,Q,Problem.b,Problem.aux.c,obj_const_term,...
                                          Problem.aux.lo,Problem.aux.hi,tol,max_IP_iter,printlevel,la_mode,fid);
        time = solution_struct.pre_time + solution_struct.runtime;
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
        fprintf(fileID,['Name = %s & IP iters = %d & & Krylov iters = %d & max_nnzL = %.2e & ' ...
            'Time = %.2e & opt = %s  \n'],Problem.name, solution_struct.IP_iter, solution_struct.Krylov_its,...
                                         solution_struct.max_nnzL, time, string(solution_struct.opt == 1)); 
        solution_statistics.objective_values(it_counter) = solution_struct.obj_val;
        solution_statistics.status_vec(it_counter) = solution_struct.opt;
        solution_statistics.problem_attempts = solution_statistics.problem_attempts + 1;
        solution_statistics.problem_names(it_counter) = string(Problem.name);
        solution_statistics.success_rate = ...
            solution_statistics.problems_converged/solution_statistics.problem_attempts;
    end
    solution_statistics.solution_struct= solution_struct;
    fprintf(fileID,['The total IPM iterates were: %d (with %d total Krylov iterates) and the total time was %d.',...
                    ' Problems converged: %d.\n'],solution_statistics.total_iters,solution_statistics.total_Krylov_iters,...
                    solution_statistics.total_time,solution_statistics.problems_converged);
    fclose(fileID);
end

