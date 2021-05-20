function [parsvec] = realpfunc(m)
%REALPFUNC Summary of this function goes here
%   Detailed explanation goes here

% parsvec(1,1) = m.I_b;   % Central moment of inertia of the chassis body about axis 3 [kg*m^2]
% parsvec(2,1) = m.I_p; % Central moment of inertia of the platform body about axis 3" [kg*m^2]
% parsvec(3,1) = m.I_a; % Axial moment of inertia of one wheel [kg*m^2]
% parsvec(4,1) = m.I_t; % Twisting moment of inertia of one wheel [kg*m^2]
% 
% parsvec(5,1) = m.l_1; % Pivot offset relative to the wheels axis [m]
% parsvec(6,1) = m.l_2; % One half of the wheels separation [m]
% 
% parsvec(7,1) = m.m_b; % Mass of the chassis base [kg]
% parsvec(8,1) = m.m_w; % Mass of one wheel [kg]
% parsvec(9,1) = m.m_p; % Mass of the platform [kg]
% 
% parsvec(10,1) = m.x_G; % x coord of the c.o.m. of the chassis body in the chassis frame [m]
% parsvec(11,1) = m.y_G; % y coord of the c.o.m. of the chassis body in the chassis frame [m]
% 
% parsvec(12,1) = m.x_F; % x coord of the c.o.m. of the platform body in the platform frame [m]
% parsvec(13,1) = m.y_F; % y coord of the c.o.m. of the platform body in the platform frame [m]
% 
% parsvec(14,1) = m.r;   % Wheel radius [m] 
% 
% parsvec(15,1) = m.b_r;   % Viscous friction coefficient for right wheel [kg*m^2*s^-1] 
% parsvec(16,1) = m.b_l;   % Viscous friction coefficient for left wheel [kg*m^2*s^-1]  
% parsvec(17,1) = m.b_p;   % Viscous friction coefficient for pivot joint [kg*m^2*s^-1] 
end

