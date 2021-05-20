%% This is the main script to call the oders and perform de model estimation

% Preset which parameters are fixed and which will be estimated
FixorNot(1,1) = false;   % I [kg*m^2]     % Paramater we want to estimate
FixorNot(2,1) = false;   % b [kg*m^3/s^2] % Parameter we want to estimate

%% Set the directories for saving the results

DirectoryName1 = 'Wheel_Motor/';
% DirectoryName1 = 'Pivot_Motor/';

% DirectoryName2 = '1.01 - Cte_Torque_SNoise [0.01]/';
DirectoryName2 = 'TEST/';

%% Set the name of the sampling plots

CN1_SP = '  Cte_T_1.01_0_5s_SNoise.png'; % Name the image

DirName_SP = '0.5s'; % Save the sample plots inside the desired folder

FullDirNameSP = strcat('Exp/',DirectoryName1,DirectoryName2,'Sample_Plots/',DirName_SP);

[MKdirStatusSP, msgSP] = mkdir(FullDirNameSP);

%% Set algorithm search options
NLGreySearchConfig = 0; % With this option the user decides the set of options concerning the Search Method and the algorithm
                        % (NLGreySearchConfig == 0)> Default setttings
                        % SearchMethod is lsqnonlin and the Algorithm is rust-region-reflective
                        
                        % (NLGreySearchConfig == 1)> Settings will be set
                        % to work with fmincon usinf the interior-point
                        % algorithm
                        
%% Set up the saving directory for the excel Tables of results

ExcelFileName = 'Ib';

DirectoryName3 = 'I_b 0.5s';

FullDirNameExcel = strcat('Exp/',DirectoryName1,DirectoryName2,'Results/',DirectoryName3,'/',ExcelFileName);
PartialDirNameExcel = strcat('Exp/',DirectoryName1,DirectoryName2,'Results/',DirectoryName3);
[MKdirStatusExcel,msgExcel] = mkdir(PartialDirNameExcel);

%% Runing the submain script
submain_nlgreyest