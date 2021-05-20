function [parsvec] = realpfunc(m)
%REALPFUNC Summary of this function goes here
%   Detailed explanation goes here
parsvec(1,1) = m.I;   % Moment of inertia of motorshaft [kg*m^2]
parsvec(2,1) = m.b;   % Viscous friction coefficient [kg*m^2*s^-1]

end

