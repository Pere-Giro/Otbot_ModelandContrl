%% Seccio 1
CN=' Experiment4_Amplified2 w_(5,1,0).png'; % Name the image

switch KinEnFlag
    case 'NO'
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
    case 'YES'
        for i=1:25
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
                names='13 - Total Kinetic Energy';
            elseif i==14
                names='14 - Rotational Kinetic Energy of Otbot';
            elseif i==15
                names='15 - Translational Kinetic Energy of Otbot';
            elseif i==16
                names='16 - Rotational Kinetic Energy of the Chassis';
            elseif i==17
                names='17 - Translational Kinetic Energy of the Chassis';         
            elseif i==18
                names='18 - Chassis Body Rotational Energy';
            elseif i==19
                names='19 - Right Wheel Rotational Energy';
            elseif i==20
                names='20 - Left Wheel Rotational Energy';
            elseif i==21
                names='21 - Platform Rotational Energy';
            elseif i==22
                names='22 - Chassis Body Translational Energy';
            elseif i==23
                names='23 - Right Wheel Translational Energy';
            elseif i==24
                names='24 - Left Wheel Rotational Energy';
            elseif i==25
                names='25 - Platform Translational Energy';
            end
            saveas(figure(i),strcat(names,CN));
        end
end

        





