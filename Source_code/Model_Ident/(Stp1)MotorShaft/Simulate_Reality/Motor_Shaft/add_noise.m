function out_data = add_noise(Noise_Flag, y_out_data1, SENSOR)
%ADD_NOISE_IMU Summary of this function goes here
%   This function adds gaussin noise to the samples data

mu = 0; % mean
sdv = SENSOR.Stdev;

[f1,c1] = size(y_out_data1);
uy = rand(f1,c1);
ny = zeros(f1,c1);

if strcmp(Noise_Flag,'NO')
    out_data = y_out_data1;
    
elseif strcmp(Noise_Flag,'YES')
    ny(:,1) = icdf('Normal',uy(:,1),mu,sdv);
    out_data = y_out_data1 + ny;

else
    disp('WARNING: This Noise_Flag does not exist, data will be sampled with no noise')
    out_data = y_out_data1;
end

end

