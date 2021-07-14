%clear;
%clc; 

% Include Code and Data files.
curr_path = pwd;
addpath(genpath(curr_path)); 

%problem_set = {"Netlib","Maros-Meszaros"};
fid = 1;   

fprintf('Should the default parameters be used?\n');
default = input('Type 1 for default parameters, or anything else to manually include them.\n');
if (default == 1)
    tol = 1e-6;                                             % Tolerance used.
    max_IP_iter = 200;                                      % Maximum number of IP iterations.
    printlevel = 3;                                         % Printing choice (see IP-PMM documentation).
    la_mode = "inexact";
    problem_set = "Maros-Meszaros";
    %problem_set = "Netlib";
    %problem_set = "SuiteSparse";
 %  problem_set = "Pearson_PDE_Optimization";
else
    fprintf('Choose a value for the allowed error tolerance.\n');
    while(true)
        tol = input('Type a double value in the form 1e-k, where k must be in [2,12].\n');
        if (isinf(tol) || isnan(tol) || ~isa(tol,'double') || tol > 1e-2 || tol < 1e-12)
            fprintf('Incorrect input argument.\n');
        else
            break;
        end
    end
    fprintf('Choose the maximum number of PMM iterations.\n');
    while(true)
        max_IP_iter = input('Type an integer value in k between [50,300].\n');
        if (isinf(max_IP_iter) || isnan(max_IP_iter) || floor(max_IP_iter)~= max_IP_iter || ...
            max_IP_iter > 300 || max_IP_iter < 50)
            fprintf('Incorrect input argument.\n');
        else
            break;
        end
    end
    fprintf('Choose the printlevel.\n');
    fprintf('                         0: no printing\n');
    fprintf('                         1: print IP-PMM iterations\n');
    fprintf('                         2: print step-length and regularization parameters\n');
    fprintf('                         3: also print Krylov iterations\n');
    while(true)
        printlevel = input('Type an integer value in k between [1,4].\n');
        if (isinf(printlevel) || isnan(printlevel) || ...
            floor(printlevel)~= printlevel || printlevel > 4 || printlevel < 1)
            fprintf('Incorrect input argument.\n');
        else
            break;
        end
    end
    fprintf('Choose the linear algebra set-up.\n');
    fprintf('                         "exact":  use factorization\n');
    fprintf('                         "inexact: use iterative linear algebra"\n');
    while(true)
        la_mode = input('Type either "exact" or "inexact" **with** the quote marks.\n');
        if (la_mode ~= "exact" && la_mode ~= "inexact")
            fprintf('Incorrect linear algebra mode argument.\n');
        else
            break;
        end
    end
    fprintf('Choose the problem set.\n');
    fprintf('                         "Netlib":  Netlib collection of LP problems\n');
    fprintf('                         "Maros-Meszaros": collection of convex QP problems"\n');
    fprintf('                         "SuiteSparse": collection of large-scale LP problems"\n');
    fprintf('                         "Pearson_PDE_Optimization": construction of PDE optimization problems"\n');
    while(true)
        problem_set = input('Type one of the previous sets **with** the quote marks.\n');
        if (problem_set ~= "Netlib" && problem_set ~= "Maros-Meszaros"&& ...
            problem_set ~= "SuiteSparse" && problem_set ~= "Pearson_PDE_Optimization")
            fprintf('Incorrect problem collection argument.\n');
        else
            break;
        end
    end
end



if (problem_set == "Netlib")
    set = 1:24;                                          % A row vector of integers in [1,96]
    type = "standard";                                  % Use "presolved" for the presolved Netlib test set.
    solution_statistics = Netlib_library(set,type,tol,max_IP_iter,printlevel,la_mode,fid);
elseif (problem_set == "Maros-Meszaros")
    set = 12:15;                                        % A row vector of integers in [1,122]
    type = "standard";                                  % Use "CONT" (1:5) for the CONT problem test set, or "standard" for the standard set.
    solution_statistics = Maros_Meszaros_library(set,type,tol,max_IP_iter,printlevel,la_mode,fid);
elseif (problem_set == "SuiteSparse")
    set = 10;
    solution_statistics = SuiteSparse_library(set,tol,max_IP_iter,printlevel,la_mode,fid);
elseif (problem_set == "Pearson_PDE_Optimization")
    % User specification of the problem to be solved.
    disp('Problem:');
    disp('         1 - Poisson Control: L^2-regularizer.');
    disp('         2 - Poisson Control: L^2-regularizer and bounded control.');
    disp('         3 - Convection Diffusion: H^1-regularizer and bounded state and control.');
    disp('         4 - Poisson Control: L^1 + L^2-regularizer and bounded control.');
    disp('         5 - Convection Diffusion: L^1 + L^2-regularizer and bounded state and control.');
    problem_choice = input('Type 1 to 5 or anything else to exit.\n');
    solution_statistics = Pearson_PDE_Test_Generator(problem_choice,tol,max_IP_iter,printlevel,la_mode,fid);
end

