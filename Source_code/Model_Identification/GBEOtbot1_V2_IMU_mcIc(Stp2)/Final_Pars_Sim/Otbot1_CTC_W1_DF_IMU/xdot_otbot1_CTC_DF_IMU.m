function [ xdot ] = xdot_otbot1_CTC_DF_IMU(t, xs, m, sm, FrictFlag, DistFlag, ulist, SampInst, Tsaction)
% XDOT_HEXAPOLE Summary of this function goes here
%   differential equation (non-linearized), that modelises our hexapole
%   system

% x = xs(1);
% y = xs(2);
alpha = xs(3);
% varphi_r = xs(4);
% varphi_l = xs(5);
varphi_p = xs(6);

% x_dot = xs(7);
% y_dot = xs(8);
alpha_dot = xs(9);
varphi_dot_r = xs(10);
varphi_dot_l = xs(11);
varphi_dot_p = xs(12);

% Setting fricctions
switch FrictFlag
    case 'YES'
        tau_friction(1:3,1) = zeros(3,1);
        tau_friction(4,1) = -m.b_r*varphi_dot_r; % tau_friction_right
        tau_friction(5,1) = -m.b_l*varphi_dot_l; % tau_friction_left
        tau_friction(6,1) = -m.b_p*varphi_dot_p; % tau_friction_pivot
    case 'NO'
        tau_friction = zeros(6,1);
    otherwise
        disp('Warning: This FrictFlag does not exist. Sistem will be simulated without friction')
        tau_friction = zeros(6,1);
end

% Setting Distrurbances
switch DistFlag
    case 'YES'
        D_vec = dist_funct(t);
    case 'NO'
        D_vec = zeros(6,1);
    otherwise
        D_vec = zeros(6,1);
end
     
Ms = sm.Mmatrix(alpha,varphi_p);
Cs = sm.Cmatrix(alpha,alpha_dot,varphi_p,varphi_dot_p);
Js = sm.Jmatrix(alpha,varphi_p);
Es = sm.Ematrix;
Jdots = sm.Jdotmatrix(alpha,alpha_dot,varphi_p,varphi_dot_p);

qdot = xs(7:12);

u_f = u_function(t, ulist, SampInst, Tsaction);
u_zoh = u_f;
qdotdot = [eye(6),zeros(6,3)]*pinv([Ms,Js.'; Js, zeros(3,3)])*[Es*u_zoh + tau_friction + D_vec - Cs*qdot; -Jdots*qdot];

xdot=[qdot; qdotdot];

end





