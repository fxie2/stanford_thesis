function cPlus = cPlus(P,tSch)
% Routine for querying c+ value for a given polynomial order
% ---
% P  - polynomial order
% tSch - time scheme - 1:RK44, 2:RK45
% cPlus - c+ value
% -------------------------------------------------------------------------

switch P
    case 2
        switch tSch
            case 1
                cPlus = 0.183;
            case 2
                cPlus = 0.206;
        end
    case 3
        switch tSch
            case 1
                cPlus = 3.6e-3;
            case 2
                cPlus = 3.8e-3;
        end
    case 4
        switch tSch
            case 1
                cPlus = 4.67e-5;
            case 2
                cPlus = 4.67e-5;
        end
    case 5
        switch tSch
            case 1
                cPlus = 4.28e-7;
            case 2
                cPlus = 4.28e-7;
        end
    otherwise
         error('cPlus:Plimit', 'Routine only provides for 2<=P<=5')
end
