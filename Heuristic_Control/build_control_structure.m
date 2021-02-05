function ctr_struct = build_control_structure()
% ==================================================================================================================== %
% build_control_structure: Store all relevant information for the heuristics controling (monitoring) IP-PMM.
% --------------------------------------------------------------------------------------------------------------------- %
% ctr_struct = build_control_structure() returns a MATLAB struct that holds the relevant information 
%                                         for controlling and monitoring the progress and state of IP-PMM.
% 
% Author: Spyridon Pougkakiotis.
% _____________________________________________________________________________________________________________________ %
    % ================================================================================================================= %
    % Store all the relevant information required for control heuristics and monitoring
    % ----------------------------------------------------------------------------------------------------------------- %
    ctr_struct = struct();
    ctr_struct.IP_iter = 0;                     % Iteration counter for IP-PMM.
    ctr_struct.opt = false;                     % Optimality indicator of the current solution. 
    ctr_struct.retry_p = 0;                     % Counter on attempts to solve the predictor.
    ctr_struct.retry_c = 0;                     % Counter on attempts to solve the corrector.
    ctr_struct.max_tries = 10;                  % Maximum number of tries: ill-conditioning message.
    ctr_struct.no_dual_update = 0;              % Primal infeasibility detection counter.
    ctr_struct.no_primal_update = 0;            % Dual infeasibility detection counter.
    ctr_struct.warn_no_est_update = 5;          % Check infeasibility if the estimates have not 
                                                % been updated for 5 iterations.
    % ________________________________________________________________________________________________________________ %
end
% ******************************************************************************************************************** %
% END OF FILE.
% ******************************************************************************************************************** %
