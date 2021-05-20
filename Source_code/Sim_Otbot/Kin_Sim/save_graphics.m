%% Seccio 1
CN=' Simulació 6 (1, 1, 1).png'; % Name the image

for i=1:6
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
    end
    saveas(figure(i),strcat(names,CN));
end





