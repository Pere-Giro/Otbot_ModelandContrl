%% Names and flags

CN=' InfinityShapeDist3'; % Name the image (Warning: avoid dots in the name if EPS save)

Flag_EPS = 'NO'; % Flag to choose if the images will be saves in eps or png)
                % (Flag_EPS == 'YES' > Images saved in EPS)
                % (Flag_EPS == 'NO'  > Images saved in png)
                
%% Saving figure

switch Flag_EPS
    case 'YES'
        saveas(figure(100),strcat(CN,'.eps'),'epsc');
    case 'NO'
        saveas(figure(100),strcat(CN,'.png'));
    otherwise
        disp('Error this Flag_EPS does not exist')
end