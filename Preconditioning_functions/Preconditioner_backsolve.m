function [w] = Preconditioner_backsolve(L_M,P,Pinv,u,D_M)
% ======================================================================== %
% [w] = Preconditioner_backsolve(L_M,P,Pinv,u,D_M): computes the backsolve
% of a Cholesky or an LDL^T decomposition.
% ------------------------------------------------------------------------ %
    if (nargin < 5 || isempty(D_M))
        D_M = [];
        Cholesky = true;
    else
        Cholesky = false;
    end
    w = u(P,1);
    if (Cholesky)
        w = L_M'\(L_M\w);
    else
        w = L_M'\(D_M\(L_M\w));
    end
    w = w(Pinv,1);
end
% ________________________________________________________________________ %
