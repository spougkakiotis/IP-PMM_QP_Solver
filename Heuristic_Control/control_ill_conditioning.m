function [ctr_struct,delta,rho,reg_limit] = control_ill_conditioning(ctr_struct,rho,delta,reg_limit,tol,predictor)
% ================================================================================================================ %
% control_ill_conditioning(ctr_struct,rho,delta,reg_limit,tol,predictor): 
% ---------------------------------------------------------------------------------------------------------------- %
% This function checks whether ill-conditioning has been detected for several iterations in a row.
% If this has happened less times than the maximum number of tries, we increase the regularization 
% parameters, as well as the regularization threshold. 
% If this occurs many times, we set the optimum status to 4, which indicates numerical inaccuracy.
% ________________________________________________________________________________________________________________ %   
    if (predictor) 
        retry = "retry_p";
    else
        retry = "retry_c";
    end
    % Checking if the matrix is too ill-conditioned. Mitigate it.
    if (ctr_struct.(retry) < ctr_struct.max_tries)
        fprintf('The system is re-solved, due to bad conditioning  of predictor system.\n')
        delta = delta*100;  rho = rho*100;
        ctr_struct.IP_iter = ctr_struct.IP_iter -1;
        ctr_struct.(retry) = ctr_struct.(retry) + 1;
        reg_limit = min(reg_limit*10,tol);
    else
        fprintf('The system matrix is too ill-conditioned.\n');
        ctr_struct.opt = 4;         % 4 signifies failure due to numerical instability.
    end
end