function PMM_struct = build_PMM_structure(zeta,lambda,delta,rho,reg_limit,mu_rate)
% ==================================================================================================================== %
% build_PMM_structure: Store all relevant information for the heuristics controling PMM parameters
% --------------------------------------------------------------------------------------------------------------------- %
% PMM_struct = build_PMM_structure(zeta,lambda,delta,rho,reg_limit,mu_rate)
%                                         returns a MATLAB struct that holds the relevant information 
%                                         for controlling and monitoring the PMM parameters.
% 
% Author: Spyridon Pougkakiotis.
% _____________________________________________________________________________________________________________________ %
    % ================================================================================================================= %
    % Store all the relevant information required for control heuristics and monitoring
    % ----------------------------------------------------------------------------------------------------------------- %
    PMM_struct = struct();
    PMM_struct.zeta = zeta;                     % Primal solution estimate.
    PMM_struct.lambda = lambda;                 % Dual multiplier estimate.
    PMM_struct.delta = delta;                   % Dual penalty (inverse) parameter.
    PMM_struct.rho = rho;                       % Primal penalty (inverse) parameter.
    PMM_struct.mu_rate = mu_rate;               % Rate of change of barrier parameter.
    PMM_struct.reg_limit = reg_limit;           % Limit for the (inverse) penalty parameters.
    % ________________________________________________________________________________________________________________ %
end
% ******************************************************************************************************************** %
% END OF FILE.
% ******************************************************************************************************************** %
