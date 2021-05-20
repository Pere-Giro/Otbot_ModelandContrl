%% This script is made to save the excel results tables

% Compute the % of deviation

percent = TSinputDev;
percent = num2str(percent,'%.0f');

fullsavename = strcat(FullDirNameExcel,percent,'%','.xlsx');

% witing the table
writetable(ResultsTable,fullsavename)
