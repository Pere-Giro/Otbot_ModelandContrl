function [X,Y,Alpha,X_dot,Y_dot,Alpha_dot,X_dotdot, Y_dotdot, Alpha_dotdot] = TD_circle_of_time_XYA(t)
%CIRCLE_OF_TIME Summary of this function goes here
%   This function is designed to compute circular trajectori for the Otbot
%   simulations

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

%-------- Experimental trajectory 3 (Infinity Shape) --------%

t1 = 4;
t2 = (3*pi)/2 + 4;
t3 = 10;
t4 = (9*pi^2-36*pi)/(8-6*pi)+10;
t5 = (12*pi-48)/(4-3*pi) + t4;

if t <= t1
    Xo = 0;
    Yo = 0;
    Vox = 0;
    Voy = 0;
    ax = 0.25;
    ay = 0;
    
    X = Xo + Vox*t + 0.5*ax*t^2;
    Y = Yo + Voy*t + 0.5*ay*t^2;
    
    X_dot = Vox + ax*t;
    Y_dot = Voy + ay*t;
    
    X_dotdot = ax; 
    Y_dotdot = ay;
    
    Alpha = 0;
    Alpha_dot = 0;
    Alpha_dotdot = 0;
    
elseif t>t1 && t<=t2
    % Xo = 2;
    % Yo = 0;
    % Vox = 1;
    % Voy = 0;
    thetao = -pi/2;
    omega = 1; % [rad/s]
    xc1 = 2;
    yc1 = 1;
    R = 1;

    X = R*cos(omega*(t-t1) + thetao) + xc1; 
    Y = R*sin(omega*(t-t1) + thetao) + yc1;

    X_dot =  - R*omega*sin(thetao + omega*(t - t1));
    Y_dot =    R*omega*cos(thetao + omega*(t - t1));

    X_dotdot =   - R*omega^2*cos(thetao + omega*(t - t1));
    Y_dotdot =   - R*omega^2*sin(thetao + omega*(t - t1));

    Alpha =  omega*(t-t1);
    Alpha_dot =  omega;
    Alpha_dotdot = 0;
    
elseif t>t2 && t<=t3
    Xo = 1;
    Yo = 1;
    Vox = 0;
    Voy = -1;
    ax = 0;
    ay = 4*(8-3*pi)/(12-3*pi)^2;
    
    X = Xo + Vox*(t-t2) + 0.5*ax*(t-t2)^2;
    Y = Yo + Voy*(t-t2) + 0.5*ay*(t-t2)^2;
    
    X_dot = Vox + ax*(t-t2);
    Y_dot = Voy + ay*(t-t2);
    
    X_dotdot = ax; 
    Y_dotdot = ay;
    
    Alpha = 3*pi/2;
    Alpha_dot = 0;
    Alpha_dotdot = 0;
    
elseif t>t3 && t<=t4
    % Xo = 1;
    % Yo = -1;
    % Vox = 0;
    % Voy = (-56+9*pi)/(16-3*pi);
    thetao = 0;
    omega = (4-3*pi)/(12-3*pi); % [rad/s]
    xc2 = 0;
    yc2 = -1;
    R = 1;

    X = R*cos(omega*(t-t3) + thetao) + xc2; 
    Y = R*sin(omega*(t-t3) + thetao) + yc2;

    X_dot =  - R*omega*sin(thetao + omega*(t - t3));
    Y_dot =    R*omega*cos(thetao + omega*(t - t3));

    X_dotdot =   - R*omega^2*cos(thetao + omega*(t - t3));
    Y_dotdot =   - R*omega^2*sin(thetao + omega*(t - t3));

    Alpha =  omega*(t-t3) + 3*pi/2;
    Alpha_dot =  omega;
    Alpha_dotdot = 0;
    
elseif t>t4 && t<=t5
    Xo = 0;
    Yo = 0;
    Vox = -((4-3*pi)/(12-3*pi));
    Voy = 0;
    Xf = 2;
    ax = 2*(Xf - Xo - Vox*(t5-t4))/(t5-t4)^2;
    ay = 0;
    
    X = Xo + Vox*(t-t4) + 0.5*ax*(t-t4)^2;
    Y = Yo + Voy*(t-t4) + 0.5*ay*(t-t4)^2;
    
    X_dot = Vox + ax*(t-t4);
    Y_dot = Voy + ay*(t-t4);
    
    X_dotdot = ax; 
    Y_dotdot = ay;
    
    Alpha = 0;
    Alpha_dot = 0;
    Alpha_dotdot = 0;
else
    X = 2;
    Y = 0;
    
    X_dot = 0;
    Y_dot = 0;
    
    X_dotdot = 0; 
    Y_dotdot = 0;
    
    Alpha = 0;
    Alpha_dot = 0;
    Alpha_dotdot = 0;
end

%--------------------------------------------------%


end
