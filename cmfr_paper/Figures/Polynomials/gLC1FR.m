function theta = gLC1FR(P,J,c1,i)

% Initialize
s = P + 3;
Cf = zeros(s+1,s+1);
rhs = zeros(s+1,1);

% Coefficients for equalities
for j = 1: P
    for k = 0: s
        Cf(j,k+1)  = (1-(-1)^(k+j))/(k+j) + (c1/J^2)*k*(1-(-1)^(k+j-2))*(j-1)/(k+j-2+eps);
    end
    rhs(j) = 0;
end

% Coefficients for constraints
for k = 0: s
    Cf(P+1,k+1)  = (-1)^k;
    Cf(P+2,k+1)  = k*(-1)^(k-1);
    Cf(P+3,k+1)  = 1;
    Cf(P+4,k+1)  = k;
    
    switch i
        case 0
            rhs(P+1) = 1;
            rhs(P+2) = 0;
            rhs(P+3) = 0;
            rhs(P+4) = 0;
        case 1
            rhs(P+1) = 0;
            rhs(P+2) = J;
            rhs(P+3) = 0;
            rhs(P+4) = 0;
    end
end


% Solve for theta
theta  = Cf \ rhs;

end