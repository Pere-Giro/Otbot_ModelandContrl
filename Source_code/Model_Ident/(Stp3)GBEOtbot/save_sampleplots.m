%% Seccio 1
NameExistence = exist('CN1_SP', 'var');

if NameExistence == 0
    CN = '  sampled_data_plots.png'; % Name the image
else
    CN = CN1_SP;
end

Direxistance = exist('FullDirNameSP','var');

if Direxistance == 0
    for i=1:5
        if i==1
            names='1 - alphadot';
        elseif i==2
            names='2 - x_ddot';
        elseif i==3
            names='3 - y_ddot';
        elseif i==4
            names='4 - taus wheels';
        elseif i==5
            names='5 - tau pivot';
        end
        saveas(figure(i),strcat(names,CN));
    end
else
    for i=1:5
        if i==1
            names='1 - alphadot';
        elseif i==2
            names='2 - x_ddot';
        elseif i==3
            names='3 - y_ddot';
        elseif i==4
            names='4 - taus wheels';
        elseif i==5
            names='5 - tau pivot';
        end
        saveas(figure(i),strcat(FullDirNameSP,'/',names,CN));
    end
end



