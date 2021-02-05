function IP_PMM_output(fid,pl,it,xinf,sinf,mu,alpha_x,alpha_z)
% =========================================================== %
% This function outputs the residual infeasibilities
% and other statistics of IP-PMM, to a give FID (file
% indicator).
% ----------------------------------------------------------- %
    if (pl >= 1)
        fprintf(fid,' ');
        fprintf(fid,'%4d    ', it);
        fprintf(fid,'%8.2e  ', xinf);
        fprintf(fid,'%8.2e  ', sinf);
        fprintf(fid,'%8.2e  ', mu);
    end
    if (pl >= 2)
        fprintf(fid,'  ');
        fprintf(fid,'%8.2e  ', alpha_x);
        fprintf(fid,'%8.2e  ', alpha_z);
    end
    if (pl >= 1) fprintf(fid,'\n'); end
% ___________________________________________________________ %
end
