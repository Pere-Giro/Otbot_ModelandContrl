%% This is the main script to call the oders and perform de model estimation

% Preset which parameters are fixed and which will be estimated
FixorNot(1,1) = true;   % I_b [kg*m^2] 
FixorNot(2,1) = false;  % I_p [kg*m^2]  % Paramater we want to estimate
FixorNot(3,1) = true;   % I_a [kg*m^2]
FixorNot(4,1) = true;   % I_t [kg*m^2]
FixorNot(5,1) = true;   % l_1 [m]
FixorNot(6,1) = true;   % l_2 [m]
FixorNot(7,1) = true;   % m_b [kg] 
FixorNot(8,1) = true;   % m_w [kg]
FixorNot(9,1) = false;  % m_p [kg]  % Paramater we want to estimate
FixorNot(10,1) = true;  % x_G [m]  
FixorNot(11,1) = true;  % y_G [m]  
FixorNot(12,1) = false; % x_F [m] % Paramater we want to estimate
FixorNot(13,1) = false; % y_F [m] % Paramater we want to estimate
FixorNot(14,1) = true;  % r   [m]
FixorNot(15,1) = true; % Abias [nounit] 
FixorNot(16,1) = true;  % Xbias [nounit] 
FixorNot(17,1) = true;  % Ybias [nounit] 

%% Set the max deviation values for each parameter
% We will be needing the radius of the platform
r_p = 0.45;           % [m]

IGmax.Ib = 0;         % [kg*m^2]
IGmax.Ip = 500*r_p^2; % [kg*m^2] % This value will not be used due to the fact that we will be computing the Ip deviation from the corresponding values of mp and xF yF
IGmax.Ia = 0;         % [kg*m^2]
IGmax.It = 0;         % [kg*m^2]

IGmax.l1 = 0;         % [m]
IGmax.l2 = 0;         % [m]

IGmax.mb = 0;         % [kg]
IGmax.mw = 0;         % [kg]
IGmax.mp = 500;       % [kg]

IGmax.x_G = 0;        % [m]
IGmax.y_G = 0;        % [m]

IGmax.x_F = r_p;      % [m]
IGmax.y_F = r_p;      % [m]

IGmax.r = 0;          % [m]

IGmax.Abias = 0;      % Nounit
IGmax.Xbias = 0;      % Nounit
IGmax.Ybias = 0;      % Nounit

%% Set the directories for saving the results

DirectoryName1 = 'Default(CoM0-0.5l1)/';
DirectoryName2 = 'Testing/';

%% Set the name of the sampling plots

CN1_SP = '  Cte_T(6,-10,6)_1.01_3s_3SNoise.png'; % Name the image

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

ExcelFileName = 'mpIpxFyF';

% DirectoryName3 = 'mp_Ip_xF_yF (1s)';
DirectoryName3 = 'mp_Ip_xF_yF (3s)';

FullDirNameExcel = strcat('Exp/',DirectoryName1,DirectoryName2,'Results/',DirectoryName3,'/',ExcelFileName);
PartialDirNameExcel = strcat('Exp/',DirectoryName1,DirectoryName2,'Results/',DirectoryName3);
[MKdirStatusExcel,msgExcel] = mkdir(PartialDirNameExcel);

%% Runing the submain script
submain_nlgreyest