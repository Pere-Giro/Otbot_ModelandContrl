function [u_out] = zoh_function(t,u)
%ZOH_FUNCTION Summary of this function goes here
%   Function designed to mantaint action signial during time steps

global tstep

if t < eps(0)
    tstep.u = u;
end

aux1_ufun = tstep.t + tstep.haction;

if t < aux1_ufun
    u_out = tstep.u;
else
    tstep.t = t;
    tstep.u = u;
    u_out = tstep.u;
end
end

