%% Seccio 1
CN=' Cte_T(6,-10,6)Data1.02_3s.png'; % Name the image

if strcmp(KinEnFlag,'NO') && strcmp(UvecFlag,'NO') && strcmp(ErrorFlag,'NO')
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
            names='13 - xddot';
        elseif i==14
            names='14 - yddot';
        elseif i==15
            names='15 - KinematicError';
        end
        saveas(figure(i),strcat(names,CN));
    end
        
elseif strcmp(KinEnFlag,'YES') && strcmp(UvecFlag,'NO') && strcmp(ErrorFlag,'NO')
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
            names='13 - xddot';
        elseif i==14
            names='14 - yddot';
        elseif i==15
            names='15 - Kinetic_Energy';
        elseif i==16
            names='16 - KinematicError';
        end
        saveas(figure(i),strcat(names,CN));
    end
elseif strcmp(KinEnFlag,'NO') && strcmp(UvecFlag,'YES') && strcmp(ErrorFlag,'NO')
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
            names='13 - xddot';
        elseif i==14
            names='14 - yddot';
        elseif i==15
            names='15 - taus wheels';
        elseif i==16
            names='16 - tau pivot';
        elseif i==17
            names='17 - KinematicError';
        end
        saveas(figure(i),strcat(names,CN));
    end
elseif strcmp(KinEnFlag,'YES') && strcmp(UvecFlag,'YES') && strcmp(ErrorFlag,'NO')
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
            names='13 - xddot';
        elseif i==14
            names='14 - yddot';
        elseif i==15
            names='15 - taus wheels';
        elseif i==16
            names='16 - tau pivot';
        elseif i==17
            names='17 - Kinetic_Energy';
        elseif i==18
            names='18 - KinematicError';
        end
        saveas(figure(i),strcat(names,CN));
    end
elseif strcmp(KinEnFlag,'NO') && strcmp(UvecFlag,'NO') && strcmp(ErrorFlag,'YES')
    for i=1:21
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
            names='13 - xddot';
        elseif i==14
            names='14 - yddot';
        elseif i==15
            names='15 - Error_X';
        elseif i==16
            names='16 - Error_Y';
        elseif i==17
            names='17 - Error_Alpha';
        elseif i==18
            names='18 - Error_Xdot';
        elseif i==19
            names='19 - Error_Ydot';    
        elseif i==20
            names='20 - Error_Alphadot';
        elseif i==21
            names='21 - KinematicError';
        end
        saveas(figure(i),strcat(names,CN));
    end
elseif strcmp(KinEnFlag,'YES') && strcmp(UvecFlag,'NO') && strcmp(ErrorFlag,'YES')
    for i=1:22
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
            names='13 - xddot';
        elseif i==14
            names='14 - yddot';
        elseif i==15
            names='15 - Kinetic_Energy';
        elseif i==16
            names='16 - Error_X';
        elseif i==17
            names='17 - Error_Y';
        elseif i==18
            names='18 - Error_Alpha';
        elseif i==19
            names='19 - Error_Xdot';
        elseif i==20
            names='20 - Error_Ydot';    
        elseif i==21
            names='21 - Error_Alphadot';
        elseif i==22
            names='22 - KinematicError';
        end
        saveas(figure(i),strcat(names,CN));
    end
elseif strcmp(KinEnFlag,'NO') && strcmp(UvecFlag,'YES') && strcmp(ErrorFlag,'YES')
    for i=1:23
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
            names='13 - xddot';
        elseif i==14
            names='14 - yddot';
        elseif i==15
            names='15 - taus wheels';
        elseif i==16
            names='16 - tau pivot';
        elseif i==17
            names='17 - Error_X';
        elseif i==18
            names='18 - Error_Y';
        elseif i==19
            names='19 - Error_Alpha';
        elseif i==20
            names='20 - Error_Xdot';
        elseif i==21
            names='21 - Error_Ydot';    
        elseif i==22
            names='22 - Error_Alphadot';
        elseif i==23
            names='23 - KinematicError';
        end
        saveas(figure(i),strcat(names,CN));
    end
elseif strcmp(KinEnFlag,'YES') && strcmp(UvecFlag,'YES') && strcmp(ErrorFlag,'YES')
    for i=1:24
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
            names='13 - xddot';
        elseif i==14
            names='14 - yddot';
        elseif i==15
            names='15 - taus wheels';
        elseif i==16
            names='16 - tau pivot';
        elseif i==17
            names='17 - Kinetic_Energy';
        elseif i==18
            names='18 - Error_X';
        elseif i==19
            names='19 - Error_Y';
        elseif i==20
            names='20 - Error_Alpha';
        elseif i==21
            names='21 - Error_Xdot';
        elseif i==22
            names='22 - Error_Ydot';    
        elseif i==23
            names='23 - Error_Alphadot';
        elseif i==24
            names='24 - KinematicError';
        end
        saveas(figure(i),strcat(names,CN));
    end
    
end


