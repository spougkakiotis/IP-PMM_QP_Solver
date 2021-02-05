clear;
clc; 

% Include Code and Data files.
curr_path = pwd;
addpath(genpath(curr_path)); 
    
set = 1:96;                                         % A row vector of integers in [1,96]
type = "standard";                                  % Use "presolved" for the presolved Netlib test set.
tol = 1e-6;                                         % Tolerance used.
max_IP_iter = 200;                                  % Maximum number of IP iterations.
printlevel = 1;                                     % Printing choice (see IP-PMM documentation).
solution_statistics = Netlib_library(set,type,tol,max_IP_iter,printlevel);

