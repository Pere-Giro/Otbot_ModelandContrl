function [D] = Dist_funct(t)
%DIST_FUNCT Summary of this function goes here
%   Detailed explanation goes here

if t>5 && t<6
    D = zeros(6,1);
    D(1:2,1) = [150;0];
elseif t>10 && t<11
    D = zeros(6,1);
    D(1:2,1) = [0;150];
elseif t>15 && t<16
    D = zeros(6,1);
    D(1:2,1) = [150;0];
elseif t>20 && t<21
    D = zeros(6,1);
    D(1:2,1) = [0;-150];
else
    D = zeros(6,1);
end

end

