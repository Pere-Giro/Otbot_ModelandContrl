function [D] = Dist_funct(t)
%DIST_FUNCT Summary of this function goes here
%   Detailed explanation goes here

if t>3 && t<4
    D = zeros(6,1);
    D(1:2,1) = [0;-150];
elseif t>8 && t<9
    D = zeros(6,1);
    D(1:2,1) = [200;0];
elseif t>11 && t<12
    D = zeros(6,1);
    D(1:2,1) = [-350;0];
elseif t>30 && t<31
    D = zeros(6,1);
    D(1:2,1) = [0;-150];
else
    D = zeros(6,1);
end

end

