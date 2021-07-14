function IP_PMM_header(fid,pl,la_mode)
% =================================================================================================================== %
% header function: (for output on a give file, indicated by fid)
% pl = 1: primal-dual infeasibility and mu is printed at each iteration k
% pl = 2: primal-dual infeasibility, mu, and step-lengths and penalty parameters are printed at each iteration k
% la_mode = "exact": no Krylov method employed
% la_mode = "inexact": also print number of Krylov iterations, nnz of preconditioner, direction error
% 
% ------------------------------------------------------------------------------------------------------------------- %
    if (nargin < 3 || isempty(la_mode)) la_mode = "exact"; end
    if (pl >= 1)
        fprintf(fid,' ');
        fprintf(fid,'%4s    ', 'iter');
        fprintf(fid,'%8s  ',   'pr feas');
        fprintf(fid,'%8s  ',   'dl feas');
        fprintf(fid,'%6s  ',   'mu');
    end
    if (pl >= 2)
        fprintf(fid,'  ');
        fprintf(fid,'%8s  ', 'alpha_x');
        fprintf(fid,'%8s  ', 'alpha_z');
        fprintf(fid,'%7s  ', 'delta');
        fprintf(fid,'%7s  ', 'rho');
    end
    if (la_mode == "inexact")
        fprintf(fid,'%12s   ', 'Krylov its');
        fprintf(fid,'%9s  ', 'nnz(Prec)');
    end

    if (pl >= 1)
        fprintf(fid,'\n ====    ========  ========  ========');
    end
    if (pl >= 2)
        fprintf(fid,'  ========  ========  ========  ========');
    end
    if (la_mode == "inexact")
        fprintf(fid,'  ===========  =========');
    end
    if (pl >= 1) fprintf(fid,'\n'); end
% ___________________________________________________________________________________________________________________ %
end
