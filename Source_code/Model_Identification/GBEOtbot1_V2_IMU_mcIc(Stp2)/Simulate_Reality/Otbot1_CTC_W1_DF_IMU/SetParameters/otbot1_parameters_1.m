%% Using function to_file to create the file with otbot physical parameters
m.I_b = 1.06458;   % Central moment of inertia of the chassis body about axis 3 [kg*m^2]
m.I_p = 2.22223; % Central moment of inertia of the platform body about axis 3" [kg*m^2]
m.I_a = 1.03570*1e-2; % Axial moment of inertia of one wheel [kg*m^2]
m.I_t = 5.61007*1e-3; % Twisting moment of inertia of one wheel [kg*m^2]

m.l_1 = 0.25; % Pivot offset relative to the wheels axis [m]
m.l_2 = 0.20; % One half of the wheels separation [m]

m.m_b = 105.00; % Mass of the chassis base [kg]
m.m_w = 2.0714; % Mass of one wheel [kg]
m.m_p = 21.94795; % Mass of the platform [kg]

m.x_G = -m.l_1/2; % x coord of the c.o.m. of the chassis body in the chassis frame [m]
m.y_G = 0;       % y coord of the c.o.m. of the chassis body in the chassis frame [m]

m.x_F = 0; % x coord of the c.o.m. of the platform body in the platform frame [m]
m.y_F = 0; % y coord of the c.o.m. of the platform body in the platform frame [m]

m.r = 0.10;   % Wheel radius [m] % No 

% Dynamic friction coefficient
m.b_r = 0.18; % Viscous friction coefficient for right wheel [kg*m^2*s^-1]
m.b_l = 0.18; % Viscous friction coefficient for left wheel [kg*m^2*s^-1]
m.b_p = 0.24; % Viscous friction coefficient for pivot joint [kg*m^2*s^-1]

% Parameters to Draw

% Matrix to define the body in relative frame with the pivot point at the
% origin
m.CBprel = [0, -m.l_1, -m.l_1;
          0, -m.l_2, +m.l_2;
          zeros(1,3)];
      
% Matrix to define the body in relative frame with the pivot point at the
% origin
       
m.RWBprel = [ + m.r,        + m.r,        - m.r,        - m.r ;
           + (m.l_2)/6, - (m.l_2)/6,  - (m.l_2)/6, + (m.l_2)/6;
           zeros(1,4)];
       
% Generic rotation matrix Rz

syms ang dang

Rzmat = [cos(ang), -sin(ang), 0;
         sin(ang), cos(ang),  0;
          0,          0,      1];
      
m.Rzmat = matlabFunction(Rzmat);

% Rotation matrices for the IMU readings

Rzmat22 = [cos(ang), -sin(ang);
         sin(ang), cos(ang)];
     
m.Rzmat22 = matlabFunction(Rzmat22);
     
IRzmat22 = inv(Rzmat22);

m.IRzmat22 = matlabFunction(IRzmat22);

dIRzmat22 = [-sin(ang)*dang, cos(ang)*dang;
             -cos(ang)*dang, -sin(ang)*dang];
         
m.dIRzmat22 = matlabFunction(dIRzmat22);


%% Saving m structure
Upath = userpath;
savedirsp1 = strcat(Upath,'\Model_Identification\GBEOtbot1_V2_IMU_mcIc(Stp2)\Simulate_Reality\Otbot1_CTC_W1_DF_IMU\');
save(strcat(savedirsp1,'m_struc.mat'),"m")
clearvars
close all
clc

%% Run the other script in order to update matrices

write_matrices_struc_2
