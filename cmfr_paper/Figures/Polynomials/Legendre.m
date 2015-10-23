%% Legendre polynomials
function Leg = Legendre(k,xi)
% ---
% Evaluates the k'th ordered Legendre polynomial at xi
% k - polynomial order 
% xi - location of evaluation in [-1,1]
% ---
if (k==-1)
    Leg = 0;
    return;
elseif (k==0)
    Leg = 1;
    return;
end
Leg = ((2*k-1)/k)*xi*Legendre(k-1,xi)-((k-1)/k)*Legendre(k-2,xi);
end