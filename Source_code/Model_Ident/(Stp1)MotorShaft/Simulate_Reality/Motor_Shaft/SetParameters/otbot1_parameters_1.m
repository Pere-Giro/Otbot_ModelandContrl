%% Using function to_file to create the file with otbot physical parameters

m.I = 1.03570*1e-2; % Axial moment of inertia of one wheel [kg*m^2]

% m.I = 2.22223; % Axial moment of inertia of the pivot shaft [kg*m^2]

% Dynamic friction coefficient
m.b = 0.18;     % Viscous friction coefficient on wheels [kg*m^2*s^-1]

% m.b = 0.24;     % Viscous friction coefficient on pivot [kg*m^2*s^-1]

%% Saving m structure
Upath = userpath;
savedirsp1 = strcat(Upath,'\Model_Identification\MotorShaft(Stp1)\Simulate_Reality\Motor_Shaft\');
save(strcat(savedirsp1,'m_struc.mat'),"m")
clearvars
close all
clc
