%% Seccio 1
CN=' Cte_T(6_-10_6)Data1_01_3s_3SNoise_25.png'; % Name the image

Flag_EPS = 'YES'; % Flag to choose if the images will be saves in eps or png)
                % (Flag_EPS == 'YES' > Images saved in EPS)
                % (Flag_EPS == 'NO'  > Images saved in png)

if strcmp(KinEnFlag,'NO') && strcmp(UvecFlag,'NO')
    for i=1:15
        if i==1
            names='1 - x';
        elseif i==2
            names='2 - y';
        elseif i==3
            names='3 - alpha';
        elseif i==4
            names='4 - varphi_r';
        elseif i==5
            names='5 - varphi_l';
        elseif i==6
            names='6 - varphi_p';
        elseif i==7
            names='7 - xdot';
        elseif i==8
            names='8 - ydot';
        elseif i==9
            names='9 - alphadot';
        elseif i==10
            names='10 - varphidot_r';
        elseif i==11
            names='11 - varphidot_l';
        elseif i==12
            names='12 - varphidot_p';
        elseif i==13
            names='13 - x_ddot';
        elseif i==14
            names='14 - y_ddot';
        elseif i==15
            names='15 - KinematicError';
        end
        switch Flag_EPS
            case 'YES'
                saveas(figure(i),strcat(names,CN,'.eps'),'epsc');
            case 'NO'
                saveas(figure(i),strcat(names,CN,'.png'));
            otherwise
                disp('Error this Flag_EPS does not exist')
        end
    end
        
elseif strcmp(KinEnFlag,'YES') && strcmp(UvecFlag,'NO')
    for i=1:16
        if i==1
            names='1 - x';
        elseif i==2
            names='2 - y';
        elseif i==3
            names='3 - alpha';
        elseif i==4
            names='4 - varphi_r';
        elseif i==5
            names='5 - varphi_l';
        elseif i==6
            names='6 - varphi_p';
        elseif i==7
            names='7 - xdot';
        elseif i==8
            names='8 - ydot';
        elseif i==9
            names='9 - alphadot';
        elseif i==10
            names='10 - varphidot_r';
        elseif i==11
            names='11 - varphidot_l';
        elseif i==12
            names='12 - varphidot_p';
        elseif i==13
            names='13 - x_ddot';
        elseif i==14
            names='14 - y_ddot';
        elseif i==15
            names='15 - Kinetic_Energy';
        elseif i==16
            names='16 - KinematicError';
        end
        switch Flag_EPS
            case 'YES'
                saveas(figure(i),strcat(names,CN,'.eps'),'epsc');
            case 'NO'
                saveas(figure(i),strcat(names,CN,'.png'));
            otherwise
                disp('Error this Flag_EPS does not exist')
        end
    end
elseif strcmp(KinEnFlag,'NO') && strcmp(UvecFlag,'YES') 
    for i=1:17
        if i==1
            names='1 - x';
        elseif i==2
            names='2 - y';
        elseif i==3
            names='3 - alpha';
        elseif i==4
            names='4 - varphi_r';
        elseif i==5
            names='5 - varphi_l';
        elseif i==6
            names='6 - varphi_p';
        elseif i==7
            names='7 - xdot';
        elseif i==8
            names='8 - ydot';
        elseif i==9
            names='9 - alphadot';
        elseif i==10
            names='10 - varphidot_r';
        elseif i==11
            names='11 - varphidot_l';
        elseif i==12
            names='12 - varphidot_p';
        elseif i==13
            names='13 - x_ddot';
        elseif i==14
            names='14 - y_ddot';
        elseif i==15
            names='15 - taus wheels';
        elseif i==16
            names='16 - tau pivot';
        elseif i==17
            names='17 - KinematicError';
        end
        switch Flag_EPS
            case 'YES'
                saveas(figure(i),strcat(names,CN,'.eps'),'epsc');
            case 'NO'
                saveas(figure(i),strcat(names,CN,'.png'));
            otherwise
                disp('Error this Flag_EPS does not exist')
        end
    end
elseif strcmp(KinEnFlag,'YES') && strcmp(UvecFlag,'YES') 
    for i=1:18
        if i==1
            names='1 - x';
        elseif i==2
            names='2 - y';
        elseif i==3
            names='3 - alpha';
        elseif i==4
            names='4 - varphi_r';
        elseif i==5
            names='5 - varphi_l';
        elseif i==6
            names='6 - varphi_p';
        elseif i==7
            names='7 - xdot';
        elseif i==8
            names='8 - ydot';
        elseif i==9
            names='9 - alphadot';
        elseif i==10
            names='10 - varphidot_r';
        elseif i==11
            names='11 - varphidot_l';
        elseif i==12
            names='12 - varphidot_p';
        elseif i==13
            names='13 - x_ddot';
        elseif i==14
            names='14 - y_ddot';
        elseif i==15
            names='15 - taus wheels';
        elseif i==16
            names='16 - tau pivot';
        elseif i==17
            names='17 - Kinetic_Energy';
        elseif i==18
            names='18 - KinematicError';
        end
        switch Flag_EPS
            case 'YES'
                saveas(figure(i),strcat(names,CN,'.eps'),'epsc');
            case 'NO'
                saveas(figure(i),strcat(names,CN,'.png'));
            otherwise
                disp('Error this Flag_EPS does not exist')
        end
    end
end


