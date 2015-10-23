function theta = gLC2FR(P,J,c,i)

% Arguments
c1 = c(1);
c2 = c(2);

% Initialize
s = P + 5;
Cf = zeros(s+1,s+1);
rhs = zeros(s+1,1);

% Coefficients for equalities
for j = 1: P
    for k = 0: s
        Cf(j,k+1)  = j*(1-(-1)^(k+j))/(k+j) + ...
            (c1/J^2)* (j*(j-1))* k* (1-(-1)^(k+j-2))/(k+j-2+eps) + ...
            (c2/J^4)* (j*(j-1)*(j-2))* (k*(k-1))* (1-(-1)^(k+j-4))/(k+j-4+eps);
    end
    rhs(j) = 0;
end

% Coefficients for constraints
for k = 0: s
    Cf(P+1,k+1)  = (-1)^k;
    Cf(P+2,k+1)  = k*(-1)^(k-1);
    Cf(P+3,k+1)  = k*(k-1)*(-1)^(k-2);
    Cf(P+4,k+1)  = 1;
    Cf(P+5,k+1)  = k;
    Cf(P+6,k+1)  = k*(k-1);
    
    switch i
        case 0
            rhs(P+1) = 1;
            rhs(P+2) = 0;
            rhs(P+3) = 0;
            rhs(P+4) = 0;
            rhs(P+5) = 0;
            rhs(P+6) = 0;
        case 1
            rhs(P+1) = 0;
            rhs(P+2) = J;
            rhs(P+3) = 0;
            rhs(P+4) = 0;
            rhs(P+5) = 0;
            rhs(P+6) = 0;
        case 2
            rhs(P+1) = 0;
            rhs(P+2) = 0;
            rhs(P+3) = J^2;
            rhs(P+4) = 0;
            rhs(P+5) = 0;
            rhs(P+6) = 0;
    end
end


% Solve for theta
theta  = Cf \ rhs;

end