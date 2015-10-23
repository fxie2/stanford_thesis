%% Derivative of Legendre polynomials
function Leg_p = Legendre_p(k,xi)
% ---
% Evaluates the derivative of the k'th ordered Legendre polynomial at xi
% k - polynomial order 
% xi - location of evaluation in [-1,1]
% ---
if (k==0 || k==-1)
    Leg_p = 0;
    return;
end
Leg_p = ((2*k-1)/k)*(Legendre(k-1,xi)+xi*Legendre_p(k-1,xi))-((k-1)/k)*Legendre_p(k-2,xi);
end