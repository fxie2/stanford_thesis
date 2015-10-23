function [err0,err1] = errNorm(uex, uex_d, t)
% Routine for calculating the 2 norms of error in 0th and 1st derivatives
% ---
% uex   - Exact solution
% uex_d - Derivative of exact solution
% -------------------------------------------------------------------------

global nelem xedge

err0 = 0;
err1 = 0;

for j = 1: nelem
    
    % Error norm for 0th derivative
    e0 = @(x) 0.5 * (uex(x,t) - uSamp0(j,x)).^2;
    err0 = err0 + quadgk(e0,xedge(j),xedge(j+1),'AbsTol',1e-7,'RelTol',1e-7);
    
    % Solution derivative
    e1 = @(x) 0.5 * (uex_d(x,t) - uSamp1(j,x)).^2;
    err1 = err1 + quadgk(e1,xedge(j),xedge(j+1),'AbsTol',1e-7,'RelTol',1e-7);
    
end

err0 = sqrt(err0);
err1 = sqrt(err1);

end

% Polynomial solution
function u0 = uSamp0(j,xl)

global u h xedge P xi

u0 = 0;
xik = -1 + (2/h(j))*(xl-xedge(j));
for p = 1: P+1
    u0 = u0 + u(p,j)*Lagrange(xi,p,xik);
end

end

% Polynomial solution derivative
function u1 = uSamp1(j,xl)

global u h xedge P xi

u1 = 0;
xik = -1 + (2/h(j))*(xl-xedge(j));
for p = 1: P+1
    u1 = u1 + u(p,j)*Lagrange_p(xi,p,xik);
end
u1 = u1 / (h(j)/2);

end