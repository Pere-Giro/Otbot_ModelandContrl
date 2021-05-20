%% Using function to_file to create the file with otbot physical parameters
I_b = 409645.04*1e-9;   % Central moment of inertia of the chassis body about axis 3 [kg*m^2]
I_p = 11262205.59*1e-9; % Central moment of inertia of the platform body about axis 3" [kg*m^2]
I_a = 11262205.59*1e-9; % Axial moment of inertia of one wheel [kg*m^2]
I_t = 11262205.59*1e-9; % Twisting moment of inertia of one wheel [kg*m^2]

l_1 = 0.2; % Pivot offset relative to the wheels axis [m]
l_2 = 0.1; % One half of the wheels separation [m]

m_b = 0.59895; % Mass of the chassis base [kg]
m_w = 0.41472; % Mass of one wheel [kg]
m_p = 0.59895; % Mass of the platform [kg]

x_G = 0.15; % x coord of the c.o.m. of the chassis body in the chassis frame [m]
y_G = 0.04; % y coord of the c.o.m. of the chassis body in the chassis frame [m]

x_F = 0.10; % x coord of the c.o.m. of the platform body in the platform frame [m]
y_F = 0.10; % y coord of the c.o.m. of the platform body in the platform frame [m]

r = 0.05;   % Wheel radius [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Geom_Vector=[I_b, I_p, I_a, I_t, l_1, l_2, m_b, m_w, m_p, x_G, y_G, x_F, y_F, r]';

Geom_Param=to_file('geom_data_otbot',Geom_Vector, 'C:\Users\pereg\Documents\UPC\Master\TFMPractiques\MatlabWD\Sim_Otbot\Kin_Sim\'); 


%% Reading from file and putting it into geometric parameters

% Geom_Vec=from_file('geom_data_otbot.txt');
% 
% m.I_b=Geom_Vec(1,1);
% m.I_p=Geom_Vec(2,1);
% m.I_a=Geom_Vec(3,1);
% m.I_t=Geom_Vec(4,1);
% m.l_1=Geom_Vec(5,1);
% m.l_2=Geom_Vec(6,1);
% m.m_b=Geom_Vec(7,1);
% m.m_w=Geom_Vec(8,1);
% m.m_p=Geom_Vec(9,1);
% m.x_G=Geom_Vec(10,1);
% m.y_G=Geom_Vec(11,1);
% m.x_F=Geom_Vec(12,1);
% m.y_F=Geom_Vec(13,1);
% m.r=Geom_Vec(14,1);



