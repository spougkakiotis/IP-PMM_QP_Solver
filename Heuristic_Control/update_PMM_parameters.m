function [PMM_struct,ctr_struct] = update_PMM_parameters(ctr_struct,nr_res_p,new_nr_res_p,nr_res_d,new_nr_res_d,x,y,PMM_struct)
% ================================================================================================================ %
% update_PMM_parameters(ctr_struct,nr_res_p,new_nr_res_p,nr_res_d,new_nr_res_d,x,y,PMM_struct):
% ---------------------------------------------------------------------------------------------------------------- %
% Using information on the primal and dual non-regularized residuals, this function decides on 
% whether the current PMM solutions is good-enough or not. If yes, then we update the primal
% and dual estimates respectively (depending on primal and dual infeasibility resp.), 
% and aggresively decrease the penalty parameters. If not, the estimates stay the same,
% and the penalty parameters are reduced at a slower rate.
% ________________________________________________________________________________________________________________ %
    cond = (0.95*norm(nr_res_p) > norm(new_nr_res_p)) || 0.95*norm(nr_res_d) > norm(new_nr_res_d);
    if (cond)
        PMM_struct.lambda = y;
        PMM_struct.delta = max(PMM_struct.reg_limit,PMM_struct.delta*(1-PMM_struct.mu_rate));  
    else % Slower rate of decrease, to avoid losing centrality.
        PMM_struct.delta = max(PMM_struct.reg_limit,PMM_struct.delta*(1-0.666*PMM_struct.mu_rate));            
        ctr_struct.no_dual_update = ctr_struct.no_dual_update + 1;
    end
    %cond = 0.95*norm(nr_res_d) > norm(new_nr_res_d);
    if (cond)
        PMM_struct.zeta = x;
        PMM_struct.rho = max(PMM_struct.reg_limit,PMM_struct.rho*(1-PMM_struct.mu_rate));  
    else
        PMM_struct.rho = max(PMM_struct.reg_limit,PMM_struct.rho*(1-0.666*PMM_struct.mu_rate));              
        ctr_struct.no_primal_update = ctr_struct.no_primal_update + 1;
    end
end