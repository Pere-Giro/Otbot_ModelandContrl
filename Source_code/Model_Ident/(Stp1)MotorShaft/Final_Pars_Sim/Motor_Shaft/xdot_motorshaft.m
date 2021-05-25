function [ xdot ] = xdot_motorshaft(t, xs, m, ulist, SampInst, Tsaction)
% XDOT_HEXAPOLE Summary of this function goes here
%   differential equation (non-linearized), that modelises our hexapole
%   system

% varphi = xs(1);
varphidot = xs(2);

qdot = xs(2);

u_f = u_function(t, ulist, SampInst, Tsaction);
u_zoh = u_f;
qdotdot = (u_zoh - m.b*varphidot)/m.I;

xdot=[qdot; qdotdot];

end





