function [X,Y,Alpha,X_dot,Y_dot,Alpha_dot, X_dotdot, Y_dotdot, Alpha_dotdot] = TD_polyline_of_time(t)
%POLYLINE_OF_TIME Summary of this function goes here
%   Detailed explanation goes here

% Polyline type 1 verical M letter
% if t <= 5
%     X =(1/5)*t;
%     Y = 0;
%     X_dot = (1/5);
%     Y_dot = 0;
% elseif t>5 && t<=10
%     X = -(1/5)*(t-5) + 1;
%     Y = (1/5)*(t-5);
%     X_dot = -(1/5);
%     Y_dot = (1/5);
% elseif t>10 && t<=15
%     X = (1/5)*(t-10);
%     Y = (1/5)*(t-5);
%     X_dot = (1/5);
%     Y_dot = (1/5);
% elseif t>15 && t<=20
%     X = -(1/5)*(t-15) + 1;
%     Y = 2;
%     X_dot = -(1/5);
%     Y_dot = 0;
% else
%     X = 0;
%     Y = 2;
%     X_dot = 0;
%     Y_dot = 0;
% end

% Polyline Type 2 Square signal
if t <= 5
    X = (1/5)*t;
    Y = 0;
    X_dot = 1/5;
    Y_dot = 0;
elseif t>5 && t<=10
    X = 1;
    Y = (1/5)*(t-5);
    X_dot = 0;
    Y_dot = 1/5;
elseif t>10 && t<=15
    X = (1/5)*(t-10)+1;
    Y = 1;
    X_dot = 1/5;
    Y_dot = 0;
elseif t>15 && t<=20
    X = 2;
    Y = -(1/5)*(t-15)+1;
    X_dot = 0;
    Y_dot = -1/5;
elseif t>20 && t<=25
    X = (1/5)*(t-20)+2;
    Y = 0;
    X_dot = 1/5;
    Y_dot = 0;
else
    X = 3;
    Y = 0;
    X_dot = 0;
    Y_dot = 0;
end

X_dotdot = 0;
Y_dotdot = 0;

Alpha = 0;
Alpha_dot = 0;
Alpha_dotdot = 0;

end

