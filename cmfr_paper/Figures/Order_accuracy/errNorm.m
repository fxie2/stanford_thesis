function [err0,err1] = errNorm(u,uref)
% Routine for calculating the 2 norms of error in 0th and 1st derivatives
% ---
% uex   - Exact solution
% uex_d - Derivative of exact solution
% -------------------------------------------------------------------------

global nelem xedge D P h

 err0 = 0;
 err1 = 0;
%{
err0 = norm(reshape(u,[1 numel(u)])-reshape(uref,[1 numel(uref)]));
numel(uref)
up = D*u;
upref = D*uref;
err1 = norm(reshape(up,[1 numel(up)])-reshape(upref,[1 numel(upref)]));
%}
%%{

display(sprintf('Calculating error for P = %i, h = %.3f',P,h(1)))
tic;
for j = 1: nelem
    
    % Error norm for 0th derivative
    e0 = @(x) 0.5 * (uSamp0(j,x,uref) - uSamp0(j,x,u)).^2;
    err0 = err0 + quadgk(e0,xedge(j),xedge(j+1),'AbsTol',1e-7,'RelTol',1e-7);
    
    % Solution derivative
    e1 = @(x) 0.5 * (uSamp1(j,x,uref) - uSamp1(j,x,u)).^2;
%     e1([xedge(j) xedge(j+1)])
%     e0([xedge(j) xedge(j+1)])
    err1 = err1 + quadgk(e1,xedge(j),xedge(j+1),'AbsTol',1e-7,'RelTol',1e-7);
    
end
%}
err0 = sqrt(err0);
err1 = sqrt(err1);
t = toc;

display(sprintf('Error calculated in %.2f s',t))
end

% Polynomial solution
function u0 = uSamp0(j,xl,uval)

global  h xedge P xi

u0 = 0;
xik = -1 + (2/h(j))*(xl-xedge(j));
for p = 1: P+1
    u0 = u0 + uval(p,j)*Lagrange(xi,p,xik);
end

end

% Polynomial solution derivative
function u1 = uSamp1(j,xl,uval)

global  h xedge P xi

u1 = 0;
xik = -1 + (2/h(j))*(xl-xedge(j));
for p = 1: P+1
    u1 = u1 + uval(p,j)*Lagrange_p(xi,p,xik);
end
u1 = u1 / (h(j)/2);

end