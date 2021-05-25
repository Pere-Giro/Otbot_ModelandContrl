function [X,Y,Alpha,X_dot,Y_dot,Alpha_dot,X_dotdot, Y_dotdot, Alpha_dotdot] = TD_circle_of_time_XYA(t)
%CIRCLE_OF_TIME Summary of this function goes here
%   This function is designed to compute circular trajectori for the Otbot
%   simulations

%-------- Circle with constant omega --------%
% X = sin(0.3*t);
% Y = cos(0.3*t) - 1;
% 
% X_dot = 0.3*cos(0.3*t);
% Y_dot = -0.3*sin(0.3*t);
% 
% X_dotdot = -0.3^2*sin(0.3*t);
% Y_dotdot = -0.3^2*cos(0.3*t);
% 
% Alpha = 0;
% Alpha_dot = 0;
% Alpha_dotdot = 0;
%--------------------------------------------%

%-------- Circle with angular acceleration --------%
omegadot = 0.3; % [rad/s^2]

X = sin(omegadot*t^2);
Y = cos(omegadot*t^2) - 1;

X_dot = 2*omegadot*t*cos(omegadot*t^2);
Y_dot = - 2*omegadot*t*sin(omegadot*t^2);

X_dotdot = 2*omegadot*cos(omegadot*t^2) - 4*omegadot^2*t^2*sin(omegadot*t^2);
Y_dotdot = - 2*omegadot*sin(omegadot*t^2) - 4*omegadot^2*t^2*cos(omegadot*t^2);

Alpha = 0;
Alpha_dot = 0;
Alpha_dotdot = 0;
%--------------------------------------------------%

end
