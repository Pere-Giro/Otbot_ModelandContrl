function [ u_f ] = u_function( t, ulist, SampInst, Tsaction)
%U_FUNCTION_LAGRANGEV2 Summary of this function goes here
%   This function has the task to give the correct action torques from the
%   ulist at every instant of time matching with the SampInst array

n = fix(t/Tsaction);

timestep = n*Tsaction;

idx = 1;

for i = 1:max(size(SampInst))
    residual = SampInst(i) - timestep;
    if residual<1e-10
        idx = i;
    end
end

u_f = ulist(idx,:)';

end
