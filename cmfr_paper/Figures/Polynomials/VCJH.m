%% Left and right VCJH correction functions
function VCJHv = VCJH(P,c,xi,bdry)
% ---
% Evaluates the (P+1)'th ordered 'bdry' VCJH correction polynomial at xi
% P+1 - polynomial order 
% c - VCJH scheme parameter
% xi - location of evaluation in [-1,1]
% bdry - -1 for left VCJH polynomial, 1 for right VCJH polynomial
% ---
aP = factorial(2*P)/((2^P)*(factorial(P))^2);
etaP = c*(2*P+1)*((aP*factorial(P))^2)/2;
if (bdry==-1)
    VCJHv = (((-1)^P)/2)*(Legendre(P,xi) - ...
        (etaP*Legendre(P-1,xi)+Legendre(P+1,xi))/(1+etaP));
elseif (bdry==1)
    VCJHv = (1/2)*(Legendre(P,xi) + ...
        (etaP*Legendre(P-1,xi)+Legendre(P+1,xi))/(1+etaP));
end
end
