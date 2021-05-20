%% Seccio 1
CN=' Cte_T_Data1_01_1.5s.png'; % Name the image

for i=1:3
    if i==1
        names='1 - varphi';
    elseif i==2
        names='2 - varphidot';
    elseif i==3
        names='3 - tau';
    end
    saveas(figure(i),strcat(names,CN));
end



