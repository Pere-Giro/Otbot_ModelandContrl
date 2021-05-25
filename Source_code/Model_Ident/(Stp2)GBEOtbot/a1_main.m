%% This is the main script to call the oders and perform de model estimation

% Preset which parameters are fixed and which will be estimated
FixorNot(1,1) = false;   % I_b [kg*m^2] % Paramater we want to estimate
FixorNot(2,1) = true;    % I_p [kg*m^2]  
FixorNot(3,1) = true;    % I_a [kg*m^2]
FixorNot(4,1) = true;    % I_t [kg*m^2]
FixorNot(5,1) = true;    % l_1 [m]
FixorNot(6,1) = true;    % l_2 [m]
FixorNot(7,1) = false;   % m_b [kg] % Paramater we want to estimate
FixorNot(8,1) = true;    % m_w [kg]
FixorNot(9,1) = true;    % m_p [kg]  
FixorNot(10,1) = false;  % x_G [m]  % Paramater we want to estimate
FixorNot(11,1) = false;  % y_G [m]  % Paramater we want to estimate
FixorNot(12,1) = true;   % x_F [m] 
FixorNot(13,1) = true;   % y_F [m] 
FixorNot(14,1) = true;   % r   [m]
FixorNot(15,1) = true;   % b_r [kg*m^2*s^-1] 
FixorNot(16,1) = true;   % b_l [kg*m^2*s^-1] 
FixorNot(17,1) = true;   % b_p [kg*m^2*s^-1] 

%% Set the directories for saving the results

DirectoryName1 = 'Default(CoM0-0.5l1)/';
% DirectoryName2 = '1.02 - 3xCte_T_FND_3Snoise [0.01]/';
DirectoryName2 = 'TEST/';

%% Set the name of the sampling plots

CN1_SP = '  Cte_T(6,-10,6)_1.02_3s_3SNoise.png'; % Name the image

DirName_SP = '3s'; % Save the sample plots inside the desired folder

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

ExcelFileName = 'mcIcxByB';

DirectoryName3 = 'mc_Ic_xB_yB 3s';

FullDirNameExcel = strcat('Exp/',DirectoryName1,DirectoryName2,'Results/',DirectoryName3,'/',ExcelFileName);
PartialDirNameExcel = strcat('Exp/',DirectoryName1,DirectoryName2,'Results/',DirectoryName3);
[MKdirStatusExcel,msgExcel] = mkdir(PartialDirNameExcel);

%% Runing the submain script
submain_nlgreyest