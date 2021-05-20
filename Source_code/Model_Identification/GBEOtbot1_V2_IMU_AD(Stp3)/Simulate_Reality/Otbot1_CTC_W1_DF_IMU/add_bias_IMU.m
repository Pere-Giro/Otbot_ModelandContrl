function [y_out_data] = add_bias_IMU(Biasalphadot_Flag, Biasxddot_Flag, Biasyddot_Flag, states, IMUxy, hts, ls, IMU)
%ADD_BIAS_IMU Summary of this function goes here
%   This function is build to add proportional bias to the sampled data

Abias = IMU.Abias;
Xbias = IMU.Xbias;
Ybias = IMU.Ybias;

if strcmp(Biasalphadot_Flag,'NO') && strcmp(Biasxddot_Flag,'NO') && strcmp(Biasyddot_Flag,'NO')
    y_out_data(:,1) = states(1:hts:ls,9);
    y_out_data(:,2:3) = IMUxy(1:hts:ls,1:2);
    
elseif strcmp(Biasalphadot_Flag,'NO') && strcmp(Biasxddot_Flag,'NO') && strcmp(Biasyddot_Flag,'YES')
    y_out_data(:,1) = states(1:hts:ls,9);
    y_out_data(:,2) = IMUxy(1:hts:ls,1);
    y_out_data(:,3) = IMUxy(1:hts:ls,2) + Ybias*IMUxy(1:hts:ls,2);
    
elseif strcmp(Biasalphadot_Flag,'NO') && strcmp(Biasxddot_Flag,'YES') && strcmp(Biasyddot_Flag,'NO')
    y_out_data(:,1) = states(1:hts:ls,9);
    y_out_data(:,2) = IMUxy(1:hts:ls,1) + Xbias*IMUxy(1:hts:ls,1);
    y_out_data(:,3) = IMUxy(1:hts:ls,2);
    
elseif strcmp(Biasalphadot_Flag,'NO') && strcmp(Biasxddot_Flag,'YES') && strcmp(Biasyddot_Flag,'YES')
    y_out_data(:,1) = states(1:hts:ls,9);
    y_out_data(:,2) = IMUxy(1:hts:ls,1) + Xbias*IMUxy(1:hts:ls,1);
    y_out_data(:,3) = IMUxy(1:hts:ls,2) + Ybias*IMUxy(1:hts:ls,2);
    
elseif strcmp(Biasalphadot_Flag,'YES') && strcmp(Biasxddot_Flag,'NO') && strcmp(Biasyddot_Flag,'NO')
    y_out_data(:,1) = states(1:hts:ls,9) + Abias*states(1:hts:ls,9);
    y_out_data(:,2) = IMUxy(1:hts:ls,1);
    y_out_data(:,3) = IMUxy(1:hts:ls,2);

elseif strcmp(Biasalphadot_Flag,'YES') && strcmp(Biasxddot_Flag,'NO') && strcmp(Biasyddot_Flag,'YES')
    y_out_data(:,1) = states(1:hts:ls,9) + Abias*states(1:hts:ls,9);
    y_out_data(:,2) = IMUxy(1:hts:ls,1);
    y_out_data(:,3) = IMUxy(1:hts:ls,2) + Ybias*IMUxy(1:hts:ls,2);
    
elseif strcmp(Biasalphadot_Flag,'YES') && strcmp(Biasxddot_Flag,'YES') && strcmp(Biasyddot_Flag,'NO')
    y_out_data(:,1) = states(1:hts:ls,9) + Abias*states(1:hts:ls,9);
    y_out_data(:,2) = IMUxy(1:hts:ls,1) + Xbias*IMUxy(1:hts:ls,1);
    y_out_data(:,3) = IMUxy(1:hts:ls,2);
    
elseif strcmp(Biasalphadot_Flag,'YES') && strcmp(Biasxddot_Flag,'YES') && strcmp(Biasyddot_Flag,'YES')
    y_out_data(:,1) = states(1:hts:ls,9) + Abias*states(1:hts:ls,9);
    y_out_data(:,2) = IMUxy(1:hts:ls,1) + Xbias*IMUxy(1:hts:ls,1);
    y_out_data(:,3) = IMUxy(1:hts:ls,2) + Ybias*IMUxy(1:hts:ls,2);
    
else
    disp('WARNING: This Biasalphadot_Flag, Biasxddot_Flag and Biasyddot_Flag combination does not exist, data will be sampled with no bias')
    y_out_data(:,1) = states(1:hts:ls,9);
    y_out_data(:,2:3) = IMUxy(1:hts:ls,1:2);
    
end

end

