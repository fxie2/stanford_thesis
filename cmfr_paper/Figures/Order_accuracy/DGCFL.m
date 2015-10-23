function CFL = DGCFL(P,tSch)
% Routine for querying c+ value for a given polynomial order
% ---
% P  - polynomial order
% tSch - time scheme - 1:RK44, 2:RK45
% CFL - DG CFL
% -------------------------------------------------------------------------

switch P
    case 0
        switch tSch
            case 1
                CFL = 1.393;
            case 2
                CFL = 2.221;
        end
    case 1
        switch tSch
            case 1
                CFL = 0.464;
            case 2
                CFL = 0.679;
        end
    case 2
        switch tSch
            case 1
                CFL = 0.235;
            case 2
                CFL = 0.352;
        end
    case 3
        switch tSch
            case 1
                CFL = 0.139;
            case 2
                CFL = 0.220;
        end
    case 4
        switch tSch
            case 1
                CFL = 0.100;
            case 2
                CFL = 0.152;
        end
    case 5
        switch tSch
            case 1
                CFL = 0.068;
            case 2
                CFL = 0.110;
        end
    otherwise
        switch tSch
            case 1
                list = [0.464, 0.235, 0.139, 0.100, 0.068]';
            case 2
                list = [0.679, 0.352, 0.220, 0.152, 0.110]';
        end
        f = fittype('a + b/x + c/x^2');
        c = fit((1:5)',list, f, 'Startpoint',[0, 0, list(1)]);
        CFL = c(P);
        
        % error('CFL:Plimit', 'Routine only provides for 1<=P<=5')
end
