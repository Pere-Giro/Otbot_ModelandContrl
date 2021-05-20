%% Using function to_file to create the file with otbot physical parameters
load('final_pars.mat')
m.I = ParameterFinalValues(1); % Moment of inertia of motor shaft [kg*m^2]
m.b = ParameterFinalValues(2); % Viscous friction coefficient [kg*m^2*s^-1]

%% Saving m structure
Upath = userpath;
savedirsp1 = strcat(Upath,'\Model_Identification\MotorShaft(Stp1)\Final_Pars_Sim\Motor_Shaft\');
save(strcat(savedirsp1,'m_struc.mat'),"m")
clearvars
close all
clc
