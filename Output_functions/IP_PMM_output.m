function IP_PMM_output(fid,pl,it,xinf,sinf,mu,alpha_x,alpha_z,delta,rho,la_mode,Krylov_its,nnz_P)
% =========================================================== %
% This function outputs the residual infeasibilities
% and other statistics of IP-PMM, to a give FID (file
% indicator).
% la_mode = "exact": no Krylov method employed
% la_mode = "inexact": also print number of Krylov iterations, nnz of preconditioner, direction error
% ----------------------------------------------------------- %
    if (nargin < 11 || isempty(la_mode))
        la_mode = "exact";
        Krylov_its = [];
        nnz_P = [];
    elseif (nargin == 12)
        error('IP_PMM_output: Not enough input arguments.\n');
    end
    if (pl >= 1)
        fprintf(fid,' ');
        fprintf(fid,'%3d     ', it);
        fprintf(fid,'%8.2e  ', xinf);
        fprintf(fid,'%8.2e  ', sinf);
        fprintf(fid,'%8.2e ', mu);
    end
    if (pl >= 2)
        fprintf(fid,' ');
        fprintf(fid,'%2.2e  ', alpha_x);
        fprintf(fid,'%8.2e  ', alpha_z);
        fprintf(fid,'%8.2e  ', delta);
        fprintf(fid,'%8.2e  ', rho);
    end
    if (la_mode == "inexact")
        fprintf(fid,'    %4d      ', Krylov_its);
        fprintf(fid,'%8.2e  ', nnz_P);
    end
    if (pl >= 1) fprintf(fid,'\n'); end
% ___________________________________________________________ %
end
