%% Seccio 1
CN=' Variable_input_Experiment7.png'; % Name the image

if strcmp(KinEnFlag,'NO') && strcmp(Jdot_qdotFlag,'NO') && strcmp(ErrorFlag,'NO')
    for i=1:12
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
        end
        saveas(figure(i),strcat(names,CN));
    end
        
elseif strcmp(KinEnFlag,'YES') && strcmp(Jdot_qdotFlag,'NO') && strcmp(ErrorFlag,'NO')
    for i=1:13
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
            names='13 - Kinetic_Energy';
        end
        saveas(figure(i),strcat(names,CN));
    end
elseif strcmp(KinEnFlag,'NO') && strcmp(Jdot_qdotFlag,'YES') && strcmp(ErrorFlag,'NO')
    for i=1:13
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
            names='13 - Jqdot_0';
        end
        saveas(figure(i),strcat(names,CN));
    end
elseif strcmp(KinEnFlag,'YES') && strcmp(Jdot_qdotFlag,'YES') && strcmp(ErrorFlag,'NO')
    for i=1:14
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
            names='13 - Kinetic_Energy';
        elseif i==14
            names='14 - Jqdot_0';
        end
        saveas(figure(i),strcat(names,CN));
    end
elseif strcmp(KinEnFlag,'NO') && strcmp(Jdot_qdotFlag,'NO') && strcmp(ErrorFlag,'YES')
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
            names='13 - Error_X';
        elseif i==14
            names='14 - Error_Y';
        elseif i==15
            names='15 - Error_Alpha';
        elseif i==16
            names='16 - Error_Xdot';
        elseif i==17
            names='17 - Error_Ydot';    
        elseif i==18
            names='18 - Error_Alphadot';
        end
        saveas(figure(i),strcat(names,CN));
    end
elseif strcmp(KinEnFlag,'YES') && strcmp(Jdot_qdotFlag,'NO') && strcmp(ErrorFlag,'YES')
    for i=1:19
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
            names='13 - Kinetic_Energy';
        elseif i==14
            names='14 - Error_X';
        elseif i==15
            names='15 - Error_Y';
        elseif i==16
            names='16 - Error_Alpha';
        elseif i==17
            names='17 - Error_Xdot';
        elseif i==18
            names='18 - Error_Ydot';    
        elseif i==19
            names='19 - Error_Alphadot';
        end
        saveas(figure(i),strcat(names,CN));
    end
elseif strcmp(KinEnFlag,'NO') && strcmp(Jdot_qdotFlag,'YES') && strcmp(ErrorFlag,'YES')
    for i=1:19
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
            names='13 - Error_X';
        elseif i==14
            names='14 - Error_Y';
        elseif i==15
            names='15 - Error_Alpha';
        elseif i==16
            names='16 - Error_Xdot';
        elseif i==17
            names='17 - Error_Ydot';    
        elseif i==18
            names='18 - Error_Alphadot';
        elseif i==19
            names='19 - Jqdot_0';
        end
        saveas(figure(i),strcat(names,CN));
    end
elseif strcmp(KinEnFlag,'YES') && strcmp(Jdot_qdotFlag,'YES') && strcmp(ErrorFlag,'YES')
    for i=1:20
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
            names='13 - Kinetic_Energy';
        elseif i==14
            names='14 - Error_X';
        elseif i==15
            names='15 - Error_Y';
        elseif i==16
            names='16 - Error_Alpha';
        elseif i==17
            names='17 - Error_Xdot';
        elseif i==18
            names='18 - Error_Ydot';    
        elseif i==19
            names='19 - Error_Alphadot';
        elseif i==20
            names='20 - Jqdot_0';
        end
        saveas(figure(i),strcat(names,CN));
    end
    
end


