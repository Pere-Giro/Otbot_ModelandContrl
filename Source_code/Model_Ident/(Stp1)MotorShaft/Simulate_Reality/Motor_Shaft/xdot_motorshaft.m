function [ xdot ] = xdot_motorshaft(t, xs, m, u, u_cte, u_Flag)
% XDOT_HEXAPOLE Summary of this function goes here
%   differential equation (non-linearized), that modelises our hexapole
%   system

% varphi    = xs(1);
varphidot = xs(2);

qdot = xs(2);

switch u_Flag
    case 'CTE'
        qdotdot = (u_cte - m.b*varphidot)/m.I;
        
    case 'VAR'
        u_f = u_function(t, u);
        u_zoh = zoh_function(t, u_f);
        qdotdot = (u_zoh - m.b*varphidot)/m.I;
        
end
xdot = [qdot; qdotdot];

end





