function IP_PMM_header(fid,pl)
% =================================================================================================================== %
% header function: (for output on a give file, indicated by fid)
% pl = 1: primal-dual infeasibility and mu is printed at each iteration k
% pl = 2: primal-dual infeasibility, mu, sigma, and step-lengths are printed at each iteration k
% ------------------------------------------------------------------------------------------------------------------- %
    if (pl >= 1)
        fprintf(fid,' ');
        fprintf(fid,'%4s    ', 'iter');
        fprintf(fid,'%8s  ',   'pr feas');
        fprintf(fid,'%8s  ',   'dl feas');
        fprintf(fid,'%8s  ',   'mu');
    end
    if (pl >= 2)
        fprintf(fid,'  ');
        fprintf(fid,'%8s  ', 'alpha_x');
        fprintf(fid,'%8s  ', 'alpha_z');
    end
    if (pl >= 1)
        fprintf(fid,'\n ====    ========  ========  ========');
    end
    if (pl >= 2)
        fprintf(fid,'    ========  ========');
    end
    if (pl >= 1) fprintf(fid,'\n'); end
% ___________________________________________________________________________________________________________________ %
end
