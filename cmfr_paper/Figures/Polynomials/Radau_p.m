%% Derivative of left and right Radau polynomials
function Rad_p = Radau_p(P,xi,bdry)
% ---
% Evaluates the derivative of the (P+1)'th ordered 'bdry' Radau polynomial at xi
% P+1 - polynomial order 
% xi - location of evaluation in [-1,1]
% bdry - -1 for left Radau polynomial, 1 for right Radau polynomial
% ---
k = P+1;
if (bdry==-1)
    Rad_p = (((-1)^k)/2)*(Legendre_p(k,xi)-Legendre_p(k-1,xi));
elseif (bdry==1)
    Rad_p = (1/2)*(Legendre_p(k,xi)+Legendre_p(k-1,xi));
end
end