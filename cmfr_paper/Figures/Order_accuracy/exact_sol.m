function u = exact_sol(x, t, deriv)

global nu k0 kN xbdry

% Evaluate the number of terms: kN/k0 
n = kN / k0;

% Evaluate translated x coordinates for the exact grid
xT = (xbdry(1)+ mod((x-t)-xbdry(1),xbdry(2)-xbdry(1)) );

u = 0;
switch deriv
    
    % 0th derivative
    case 0
        for j = 1: n
            u = u + exp(-nu* (j*k0)^2 *t) *(1/n)* cos(j *k0 *xT);
        end
        
    % 1st derivative
    case 1
        for j = 1: n
            u = u - exp(-nu* (j*k0)^2 *t) *(1/n)* (j*k0) * sin(j *k0 *xT);
        end
        
end
    