function xi = GLpoints(P)
% Routine for querying Gauss-Legendre points for a given polynomial order
% ---
% P  - polynomial order
% xi - Gauss Legendre points
% -------------------------------------------------------------------------

switch P
    case 0
        xi = 0;
    case 1
        xi = [-sqrt(1/3); sqrt(1/3)];
    case 2
        xi = [-sqrt(3/5); 0; sqrt(3/5)];
    case 3
        xi = [-sqrt((3+2*sqrt(6/5))/7); -sqrt((3-2*sqrt(6/5))/7); sqrt((3-2*sqrt(6/5))/7); sqrt((3+2*sqrt(6/5))/7)];
    case 4
        xi = [-(1/3)*sqrt(5+2*sqrt(10/7)); -(1/3)*sqrt(5-2*sqrt(10/7)); 0; (1/3)*sqrt(5-2*sqrt(10/7)); (1/3)*sqrt(5+2*sqrt(10/7))];
    case 5
        xi = [-0.9324695142; -0.6612093865; -0.2386191861; 0.2386191861; 0.6612093865; 0.9324695142];
    case 6
        xi = [-0.9491079123; -0.7415311856; -0.4058451514; 0; 0.4058451514; 0.7415311856; 0.9491079123];
    otherwise
         error('GLpoints:Plimit', 'Routine only provides for P<=6')
end
