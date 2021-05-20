function [ u_f ] = u_function( t, u )
%U_FUNCTION_LAGRANGEV2 Summary of this function goes here
%   Detailed explanation goes here

if t>1
    u_f = [0;0;0];
else
    u_f = u;
end

end

