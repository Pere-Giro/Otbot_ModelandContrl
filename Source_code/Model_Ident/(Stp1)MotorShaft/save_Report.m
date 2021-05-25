%% This script is made to save the Results Report given by nlgreyest

filename = 'report_file';

ReportNLGRM = nlgrm.Report; % Create report object

percent = TSinputDev;

percent = num2str(percent,'%.0f');
Tsvalue = num2str(data.Ts,'%.3f');

% fullsavename = strcat('SensitivityData/',filename,Tsvalue,'s_',percent,'.mat');
fullsavename = strcat(PartialDirNameExcel,'/',filename,Tsvalue,'s_',percent,'.mat');

save(fullsavename,"ReportNLGRM")

%% Write a txt file with the termination condition
docname = 'term_cond';

TCnlgrm = ReportNLGRM.Termination;
OpUnlgrm = ReportNLGRM.OptionsUsed;

% fulldocname = strcat('SensitivityData/',docname,Tsvalue,'s_',percent,'.txt');
fulldocname = strcat(PartialDirNameExcel,'/',docname,Tsvalue,'s_',percent,'.txt');

fileID = fopen(fulldocname,'w');
fprintf(fileID,'%s %s\r\n','Method: ',ReportNLGRM.Method);
fprintf(fileID,'%s\r\n','-------------------------------');
fprintf(fileID,'%s\r\n','Termination Report');
fprintf(fileID,'%s\r\n','-------------------------------');
fprintf(fileID,'%s %s\r\n','WhyStop: ',TCnlgrm.WhyStop);
fprintf(fileID,'%s %.0f\r\n','Iterations: ',TCnlgrm.Iterations);
fprintf(fileID,'%s %.4e\r\n','FirstOrderOptimality: ',TCnlgrm.FirstOrderOptimality);
fprintf(fileID,'%s %.0f\r\n','FcnCount: ',TCnlgrm.FcnCount);
fprintf(fileID,'%s %s\r\n','Algorithm: ',TCnlgrm.Algorithm);
fprintf(fileID,'%s\r\n','-------------------------------');
fprintf(fileID,'%s\r\n','OptionsUsed Report');
fprintf(fileID,'%s\r\n','-------------------------------');
fprintf(fileID,'%s %s\r\n','SearchMethod: ',OpUnlgrm.SearchMethod);
fprintf(fileID,'%s %.0f\r\n','MaxIterations: ',OpUnlgrm.SearchOptions.MaxIterations);
fprintf(fileID,'%s %.4e\r\n','StepTolerance: ',OpUnlgrm.SearchOptions.StepTolerance);
fprintf(fileID,'%s %.4e\r\n','FunctionTolerance: ',OpUnlgrm.SearchOptions.FunctionTolerance);
fprintf(fileID,'%s\r\n','-------------------------------');
fprintf(fileID,'%s\r\n','Time to solve Report');
fprintf(fileID,'%s\r\n','-------------------------------');
fprintf(fileID,'%s %.4f\r\n','Elapsed time [s]: ',tEnd);
fclose(fileID);

