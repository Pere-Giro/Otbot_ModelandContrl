function [ qdot ] = qdot_kin_otbot(t, xs, sm, v)
%XDOT_HEXAPOLE Summary of this function goes here
%   differential equation (non-linearized), that modelises our hexapole
%   system

x = xs(1);
y = xs(2);
alpha = xs(3);
varphi_r = xs(4);
varphi_l = xs(5);
varphi_p = xs(6);

Mfik = sm.MFIKmatrix(alpha,varphi_p);

Maux = [Mfik;
         eye(3)];
     
qdot = Maux*v; 

end





