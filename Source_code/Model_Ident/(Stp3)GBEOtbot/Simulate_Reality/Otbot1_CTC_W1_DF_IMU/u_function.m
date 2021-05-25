function [ u_f ] = u_function( t, u, m, sm )
%U_FUNCTION_LAGRANGEV2 Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%  Step u  %%%%%%%%%
% if t>1
%     u_f = [0;0;0];
% else
%     u_f = u;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%  Ramp u  %%%%%%%%%
% u_f = u + 5*t*[1;1;0];
% u_f = u + 10*t*[1;1;0];
% u_f = u + 40*t*[1;1;0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Square signal 1s %%%%%%%%%
% Tpr = 1;       % [s] Period of the square signal for right wheel
% Tpl = Tpr*5/7; % [s] Period of the square signal for left wheel
% Tpp = Tpr*11/7; % [s] Period of the square signal for pivot
% Amp = 6; % [N*m] Amplitude of the signal
% 
% 
% wave(1,1) = sin(2*pi/Tpr*t + deg2rad(90)); % right wheel
% % wave(2,1) = sin(2*pi/Tp*t + deg2rad(90)); % left wheel
% % wave(3,1) = sin(2*pi/Tp*t + deg2rad(180)); % Pivot joint
% wave(2,1) = sin(2*pi/Tpl*t + deg2rad(90)); % left wheel
% wave(3,1) = sin(2*pi/Tpp*t + deg2rad(90)); % Pivot joint
% 
% % Torque for the right wheel
% if wave(1,1)>0
%     u_f(1,1) = Amp;
% elseif wave(1,1)<0
%     u_f(1,1) = -Amp;
% else
%     u_f(1,1)= 0;
% end
% 
% % Torque for the left wheel
% if wave(2,1)>0
%     u_f(2,1) = Amp;
% elseif wave(2,1)<0
%     u_f(2,1) = -Amp;
% else
%     u_f(2,1)= 0;
% end
% 
% % Torque for the pivot joint
% if wave(3,1)>0
%     u_f(3,1) = Amp;
% elseif wave(3,1)<0
%     u_f(3,1) = -Amp;
% else
%     u_f(3,1)= 0;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Sine waves 1s %%%%%%%%%
Tpr = 1;       % [s] Period of the square signal for right wheel
Tpl = Tpr*5/7; % [s] Period of the square signal for left wheel
Tpp = Tpr*11/7; % [s] Period of the square signal for pivot
Amp = 3; % [N*m] Amplitude of the signal
Offset = 3; % Offset to apply to the wave

u_f(1,1) = Amp*sin(2*pi/Tpr*t + deg2rad(90)) + Offset; % Right wheel torque
u_f(2,1) = Amp*sin(2*pi/Tpl*t + deg2rad(90)) + Offset; % Left wheel torque
u_f(3,1) = Amp*sin(2*pi/Tpp*t + deg2rad(90)) + Offset; % Pivot joint torque

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Circular Trajectory %%%%%%%%%

% % Parameters of the circle
% Rm = 1; % [m]
% omega = 1; % [rad/s]
% %---------------------------
% 
% % x_of_t = Rm*cos(omega*t) - m.l_1*sin(omega*t);
% % y_of_t = Rm*sin(omega*t) + m.l_1*cos(omega*t);
% alpha_of_t = 0;
% 
% x_dot_of_t = -Rm*omega*sin(omega*t) - m.l_1*omega*cos(omega*t);
% y_dot_of_t = Rm*omega*cos(omega*t) - m.l_1*omega*sin(omega*t);
% alpha_dot_of_t = 0;
% 
% x_dotdot_of_t = -Rm*omega^2*cos(omega*t) + m.l_1*omega^2*sin(omega*t);
% y_dotdot_of_t = -Rm*omega^2*sin(omega*t) - m.l_1*omega^2*cos(omega*t);
% alpha_dotdot_of_t = 0;
% 
% % p = [x_of_t; y_of_t; alpha_of_t];
% pdot = [x_dot_of_t; y_dot_of_t; alpha_dot_of_t];
% pdotdot = [x_dotdot_of_t; y_dotdot_of_t; alpha_dotdot_of_t];
% 
% varphi_p_of_t = -omega*t - pi/2;
% varphi_dot_p_of_t = -omega;
% 
% % Use of Inverse Dynamics to compute u trajectory
% 
% Mbars = sm.M_barmatrix(alpha_of_t,varphi_p_of_t);
% Cbars = sm.C_barmatrix(alpha_of_t,alpha_dot_of_t,varphi_p_of_t,varphi_dot_p_of_t);
% MIIKs = sm.MIIKmatrix(alpha_of_t,varphi_p_of_t); 
% 
% u_f = pinv(MIIKs.')*(Mbars*pdotdot + Cbars*pdot);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

