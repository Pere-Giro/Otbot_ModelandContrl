function [X,Y,Alpha,X_dot,Y_dot,Alpha_dot, X_dotdot, Y_dotdot, Alpha_dotdot] = TD_polyline_of_time(t)
%POLYLINE_OF_TIME Summary of this function goes here
%   %   This function is designed to compute trajectories for the Otbot
%   simulations

% Polyline type 1 verical M letter
% if t <= 5
%     X = (1/5*0.4)*t;
%     Y = 0;
%     X_dot = (1/5*0.4);
%     Y_dot = 0;
% elseif t>5 && t<=10
%     X = -(1/5*0.4)*(t-5) + 0.4;
%     Y = (1/5*0.4)*(t-5);
%     X_dot = -(1/5*0.4);
%     Y_dot = (1/5*0.4);
% elseif t>10 && t<=15
%     X = (1/5*0.4)*(t-10);
%     Y = (1/5*0.4)*(t-5); 
%     X_dot = 0.4*(1/5);
%     Y_dot = 0.4*(1/5);
% elseif t>15 && t<=20
%     X = -(1/5*0.4)*(t-15) + 0.4;
%     Y = 0.8;
%     X_dot = -(1/5*0.4);
%     Y_dot = 0;
% else
%     X = 0;
%     Y = 0.8;
%     X_dot = 0;
%     Y_dot = 0;
% end

% Polyline Type 2 Square signal
% if t <= 5
%     X = (1/5)*t;
%     Y = 0;
%     X_dot = 1/5;
%     Y_dot = 0;
% elseif t>5 && t<=10
%     X = 1;
%     Y = (1/5)*(t-5);
%     X_dot = 0;
%     Y_dot = 1/5;
% elseif t>10 && t<=15
%     X = (1/5)*(t-10)+1;
%     Y = 1;
%     X_dot = 1/5;
%     Y_dot = 0;
% elseif t>15 && t<=20
%     X = 2;
%     Y = -(1/5)*(t-15)+1;
%     X_dot = 0;
%     Y_dot = -1/5;
% elseif t>20 && t<=25
%     X = (1/5)*(t-20)+2;
%     Y = 0;
%     X_dot = 1/5;
%     Y_dot = 0;
% else
%     X = 3;
%     Y = 0;
%     X_dot = 0;
%     Y_dot = 0;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% X_dotdot = 0;
% Y_dotdot = 0;
% 
% Alpha = 0;
% Alpha_dot = 0;
% Alpha_dotdot = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant acceleration

% Set acceleration value
% a = 0.5; % [m/s^2]
% 
% X = 0.5*a*t^2;
% Y = 0;
% 
% X_dot = a*t;
% Y_dot = 0;
% 
% X_dotdot = 0;
% Y_dotdot = 0;
% 
% Alpha = 0;
% Alpha_dot = 0;
% Alpha_dotdot = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable acceleration (Acceleration function of time)

% Set acceleration value
a = 0.5; % [m/s^2]

X = (1/6)*a*t^3;
Y = 0;

X_dot = 0.5*a*t^2;
Y_dot = 0;

X_dotdot = 0;
Y_dotdot = 0;

Alpha = 0;
Alpha_dot = 0;
Alpha_dotdot = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
