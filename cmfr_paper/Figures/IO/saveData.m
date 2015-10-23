function saveData(str, t, uex, u, ind, iter)
% Routine for storing data
% ---
% str - File tag
% t - time
% uex - exact solution
% u - numerical solution
% ind - scheme index
% ---

global x xbdry xedge h nelem xi P xex nex plotStyle

switch ind
    case 1
        u1 = u;
        save([str,num2str(t),'.mat'],'xex','uex','u1', ...
            'x','u','xbdry','xedge','h','nelem','xi','P','xex','nex','plotStyle','iter','t');
    case 2
        u2 = u;
        save([str,num2str(t),'.mat'],'u2','-append');
    case 3
        u3 = u;
        save([str,num2str(t),'.mat'],'u3','-append');
    case 4
        u4 = u;
        save([str,num2str(t),'.mat'],'u4','-append');
    case 5
        u5 = u;
        save([str,num2str(t),'_c+.mat'],'xex','uex','u5',...
            'x','u','xbdry','xedge','h','nelem','xi','P','xex','nex','plotStyle','iter','t');
end
       
    
