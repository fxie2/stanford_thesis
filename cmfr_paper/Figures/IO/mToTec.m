function mToTec(str, t, uex, ind)
% Routine for storing data in tecplot format
% ---
% str - File tag
% t - time
% uex  - Exact solution
% ind - plot index
% -------------------------------------------------------------------------

global x u xbdry xedge h nelem xi P xex nex

% Open handle
fid = fopen([str, '_t_',num2str(t),'.dat'],'w');

% Exact solution
if (ind == 1 || ind ==5)
    xlist = reshape(xex,nex*nelem,1);
    ulist = reshape(uex,nex*nelem,1);
end

% Discrete numerical solution
xhlist = reshape(x,(P+1)*nelem,1);
uhlist = reshape(u,(P+1)*nelem,1);

% Constructed numerical solution
uel = zeros(nex,nelem);
for j = 1: nelem
    for k = 1: nex
        xik = -1 + (2/h(j))*(xex(k,j)-xedge(j));
        for p = 1: P+1
            uel(k,j) = uel(k,j) + u(p,j)*Lagrange(xi,p,xik);
        end
    end
end

% Write data to file
fprintf(fid,)

fprintf(fid, myformat, magic(4));
fclose(fid);

