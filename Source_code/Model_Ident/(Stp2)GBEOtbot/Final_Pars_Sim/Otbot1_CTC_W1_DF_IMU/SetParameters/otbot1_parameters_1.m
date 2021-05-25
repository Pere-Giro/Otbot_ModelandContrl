%% Using function to_file to create the file with otbot physical parameters
load('final_pars.mat')
m.I_b = ParameterFinalValues(1);   % Central moment of inertia of the chassis body about axis 3 [kg*m^2]
m.I_p = ParameterFinalValues(2); % Central moment of inertia of the platform body about axis 3" [kg*m^2]
m.I_a = ParameterFinalValues(3); % Axial moment of inertia of one wheel [kg*m^2]
m.I_t = ParameterFinalValues(4); % Twisting moment of inertia of one wheel [kg*m^2]

m.l_1 = ParameterFinalValues(5); % Pivot offset relative to the wheels axis [m]
m.l_2 = ParameterFinalValues(6); % One half of the wheels separation [m]

m.m_b = ParameterFinalValues(7); % Mass of the chassis base [kg]
m.m_w = ParameterFinalValues(8); % Mass of one wheel [kg]
m.m_p = ParameterFinalValues(9); % Mass of the platform [kg]

m.x_G = ParameterFinalValues(10); % x coord of the c.o.m. of the chassis body in the chassis frame [m]
m.y_G = ParameterFinalValues(11);       % y coord of the c.o.m. of the chassis body in the chassis frame [m]

m.x_F = ParameterFinalValues(12); % x coord of the c.o.m. of the platform body in the platform frame [m]
m.y_F = ParameterFinalValues(13); % y coord of the c.o.m. of the platform body in the platform frame [m]

m.r = ParameterFinalValues(14);   % Wheel radius [m] % No 

% Dynamic friction coefficient
m.b_r = ParameterFinalValues(15); % Viscous friction coefficient for right wheel [kg*m^2*s^-1]
m.b_l = ParameterFinalValues(16); % Viscous friction coefficient for left wheel [kg*m^2*s^-1]
m.b_p = ParameterFinalValues(17); % Viscous friction coefficient for pivot joint [kg*m^2*s^-1]

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
savedirsp1 = strcat(Upath,'\Model_Identification\GBEOtbot1_V2_IMU_mcIc(Stp2)\Final_Pars_Sim\Otbot1_CTC_W1_DF_IMU\');
save(strcat(savedirsp1,'m_struc.mat'),"m")
clearvars
close all
clc

%% Run the other script in order to update matrices

write_matrices_struc_2
