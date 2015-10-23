function u = init_cond(x, deriv)

global k0 kN

% Evaluate the number of terms: kN/k0 
n = kN / k0;

u = zeros(size(x));
switch deriv
    
    % 0th derivative
    case 0
        for j = 1: n
            u = u + (1/n)*cos(j *k0 *x);
        end
        
    % 1st derivative
    case 1
        for j = 1: n
            u = u - (1/n)*(j*k0)*sin(j *k0 *x);
        end
        
end
    