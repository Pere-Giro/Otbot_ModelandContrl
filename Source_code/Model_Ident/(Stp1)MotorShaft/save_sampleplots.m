%% Seccio 1
NameExistence = exist('CN1_SP', 'var');

if NameExistence == 0
    CN = '  sampled_data_plots.png'; % Name the image
else
    CN = CN1_SP;
end

Direxistance = exist('FullDirNameSP','var');

if Direxistance == 0
    for i=1:2
        if i==1
            names='1 - varphidot';
        elseif i==2
            names='2 - tau action';
        end
        saveas(figure(i),strcat(names,CN));
    end
else
    for i=1:2
        if i==1
            names='1 - varphidot';
        elseif i==2
            names='2 - tau pivot';
        end
        saveas(figure(i),strcat(FullDirNameSP,'/',names,CN));
    end
end



