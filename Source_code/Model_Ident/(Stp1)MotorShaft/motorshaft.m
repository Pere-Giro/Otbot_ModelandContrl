function [xdot,yout] = motorshaft(t, xs, u, I, b, varargin)
%OTBOT_M Summary of this function goes here
%  Function containing the equations of motion of Otbot

%------------ List of state variables ------------%
% Note how some of this state variables (the ones commented out) are not 
% used to compute the equations of motion of the system.

% varphi    = xs(1);
varphidot = xs(2);

%------------ Compute the equations ------------%

qdot = xs(2);

qdotdot = (u - b*varphidot)/I;

xdot = [qdot; qdotdot];

% Output equation 
yout = [xs(2)]; % varphidot

end


