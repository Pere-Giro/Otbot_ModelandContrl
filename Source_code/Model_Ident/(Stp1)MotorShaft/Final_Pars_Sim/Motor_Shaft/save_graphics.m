%% Seccio 1
CN=' Cte_TData1_01_0_5s_SNoise_50'; % Name the image (Warning: avoid dots in the name if EPS save)

Flag_EPS = 'YES'; % Flag to choose if the images will be saves in eps or png)
                % (Flag_EPS == 'YES' > Images saved in EPS)
                % (Flag_EPS == 'NO'  > Images saved in png)

for i=1:3
    if i==1
        names='1 - varphi';
    elseif i==2
        names='2 - varphidot';
    elseif i==3
        names='3 - tau';
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

