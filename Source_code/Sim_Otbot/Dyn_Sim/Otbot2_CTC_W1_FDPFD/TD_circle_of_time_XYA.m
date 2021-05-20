function [X,Y,Alpha,X_dot,Y_dot,Alpha_dot,X_dotdot, Y_dotdot, Alpha_dotdot] = TD_circle_of_time_XYA(t)
%CIRCLE_OF_TIME Summary of this function goes here
%   This function is designed to compute circular trajectori for the Otbot
%   simulations

X = 0.4*sin(0.3*t);
Y = 0.4*cos(0.3*t) - 0.4;

X_dot = 0.4*0.3*cos(0.3*t);
Y_dot = -0.4*0.3*sin(0.3*t);

X_dotdot = -0.4*0.3^2*sin(0.3*t);
Y_dotdot = -0.4*0.3^2*cos(0.3*t);

Alpha = 0;
Alpha_dot = 0;
Alpha_dotdot = 0;

end
