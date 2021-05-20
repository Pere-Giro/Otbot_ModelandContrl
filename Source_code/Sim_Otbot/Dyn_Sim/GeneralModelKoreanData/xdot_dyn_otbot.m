function [ xdot ] = xdot_dyn_otbot(t, xs, sm, u, u_Flag)
%XDOT_HEXAPOLE Summary of this function goes here
%   differential equation (non-linearized), that modelises our hexapole
%   system

x = xs(1);
y = xs(2);
alpha = xs(3);
varphi_r = xs(4);
varphi_l = xs(5);
varphi_p = xs(6);

x_dot = xs(7);
y_dot = xs(8);
alpha_dot = xs(9);
varphi_dot_r = xs(10);
varphi_dot_l = xs(11);
varphi_dot_p = xs(12);

Ms = sm.Mmatrix(alpha,varphi_p);
Js = sm.Jmatrix(alpha,varphi_p);
Es = sm.Ematrix;
Cs = sm.Cmatrix(alpha,alpha_dot,varphi_p,varphi_dot_p);
Jdots = sm.Jdotmatrix(alpha,alpha_dot,varphi_p,varphi_dot_p);

qdot = xs(7:12);

switch u_Flag
    case 'CTE'
        qdotdot = [eye(6),zeros(6,3)]*inv([Ms,Js.'; Js, zeros(3,3)])*[Es*u - Cs*qdot; -Jdots*qdot];
    case 'VAR'
        u_f= u_function(t,u);
        qdotdot = [eye(6),zeros(6,3)]*inv([Ms,Js.'; Js, zeros(3,3)])*[Es*u_f - Cs*qdot; -Jdots*qdot];
    otherwise
        disp('WARNING THIS FLAG DOES NOT EXIST')
        disp('Simulation with CTE torques will be launched')
        qdotdot = [eye(6),zeros(6,3)]*inv([Ms,Js.'; Js, zeros(3,3)])*[Es*u - Cs*qdot; -Jdots*qdot];
end
xdot=[qdot; qdotdot];
end





