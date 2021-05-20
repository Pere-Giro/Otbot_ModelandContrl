%% This script is made to save the required data to make the sensitivity plots

filename = 'sensitivity_data';

% Setting up the Free or fixed parameters for estimation vector

FixedVector = [nlgrm.parameters(1).Fixed;   % I_b [kg*m^2] 
                nlgrm.parameters(2).Fixed;  % I_p [kg*m^2]
                nlgrm.parameters(3).Fixed;  % I_a [kg*m^2]
                nlgrm.parameters(4).Fixed;  % I_t [kg*m^2]
                nlgrm.parameters(5).Fixed;  % l_1 [m]
                nlgrm.parameters(6).Fixed;  % l_2 [m]
                nlgrm.parameters(7).Fixed;  % m_b [kg] 
                nlgrm.parameters(8).Fixed;  % m_w [kg]
                nlgrm.parameters(9).Fixed;  % m_p [kg] 
                nlgrm.parameters(10).Fixed; % x_G [m]
                nlgrm.parameters(11).Fixed; % y_G [m]
                nlgrm.parameters(12).Fixed; % x_F [m]
                nlgrm.parameters(13).Fixed; % y_F [m]
                nlgrm.parameters(14).Fixed; % r   [m]
                nlgrm.parameters(15).Fixed; % Abias [nounit]
                nlgrm.parameters(16).Fixed; % Xbias [nounit]
                nlgrm.parameters(17).Fixed];% Ybias [nounit]
            
Dev_value = ParameterInitialGuesses;

%% Saveing all the Data in a .mat file

percent = TSinputDev;
percent = num2str(percent,'%.0f');
Tsvalue = num2str(data.Ts,'%.3f');

fullsavename = strcat('SensitivityData/',filename,Tsvalue,'s_',percent,'.mat');

save(fullsavename,"Dev_value","FixedVector","RelativeError_Percentage","AbsoluteError")

