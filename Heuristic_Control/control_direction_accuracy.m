function [ctr_struct,mu,delta,rho] = control_direction_accuracy(ctr_struct,mu,mu_prev,delta,rho,predictor)
% ================================================================================================================ %
% control_direction_accuracy(ctr_struct,mu,mu_prev,predictor): 
% ---------------------------------------------------------------------------------------------------------------- %
% This function checks whether an inaccurate direction has been encountered many times in the row.
% If the accuracy is not good, it forces IP-PMM to backtrack and re-calculate it.
% If this happens many times consequtively, it forces termination with a message indicating 
% inaccuracy. 
% ________________________________________________________________________________________________________________ %   
    if (predictor) 
        retry = "retry_p";
    else
        retry = "retry_c";
    end    
    if (ctr_struct.(retry) < ctr_struct.max_tries)
        fprintf('Corrector: Dropping the direction, due to inaccuracy.\n');
        ctr_struct.IP_iter = ctr_struct.IP_iter-1;
        ctr_struct.(retry) = ctr_struct.(retry) + 1;
        mu = mu_prev;
        delta = delta*10;
        rho = rho*10;
        ctr_struct.no_primal_update = max(ctr_struct.no_primal_update-1,0);
        ctr_struct.no_dual_update = max(ctr_struct.no_dual_update-1,0);
     else
        fprintf('Not enough accuracy.\n');
        ctr_struct.opt = 5;         % 5 signifies failure due to inaccuracy.
     end
end