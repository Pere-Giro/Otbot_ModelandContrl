function [parsvec] = realpfunc(m)
%REALPFUNC Summary of this function goes here
%   this function computes the required model parameters form the m struc


parsvec(1,1) = m.I; % Moment of inertia of the motorshaft [kg*m^2]
parsvec(2,1) = m.b; % Viscous friction coefficient [kg*m^2*s^-1]

end

