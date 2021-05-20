function out_data = add_noise_IMU(Noisealphadot_Flag, Noisexddot_Flag, Noiseyddot_Flag, Nopt_Flag1, Nopt_Flag2, Nopt_Flag3, y_out_data1, IMU)
%ADD_NOISE_IMU Summary of this function goes here
%   This function adds gaussin noise to the samples data

ASR = IMU.ASR;
XSR = IMU.XSR;
YSR = IMU.YSR;

NDgy = IMU.NDgy;
NDaccX = IMU.NDaccX;
NDaccY = IMU.NDaccX;

mu = 0; % mean

% Compute st devs
if strcmp(Nopt_Flag1,'DATA') 
    sdv(1) = NDgy*sqrt(ASR);   % Noise on Alphadot
elseif strcmp(Nopt_Flag1,'VALUE')
    sdv(1) = IMU.StdevA;
else
    disp('WARNING: This Nopt_Flag1, does not exist the Stdev of the Noise for alpha_dot will be set to the IMU.StdevA value' )
    sdv(1) = IMU.StdevA;
end

if strcmp(Nopt_Flag2,'DATA') 
    sdv(2) = NDaccX*sqrt(XSR); % Noise on Xddot
elseif strcmp(Nopt_Flag2,'VALUE')
    sdv(2) = IMU.StdevX;
else
    disp('WARNING: This Nopt_Flag1, does not exist the Stdev of the Noise for alpha_dot will be set to the IMU.StdevA value' )
    sdv(2) = IMU.StdevX;
end

if strcmp(Nopt_Flag3,'DATA') 
    sdv(3) = NDaccY*sqrt(YSR); % Noise on Yddot
elseif strcmp(Nopt_Flag3,'VALUE')
    sdv(3) = IMU.StdevY;
else
    disp('WARNING: This Nopt_Flag1, does not exist the Stdev of the Noise for alpha_dot will be set to the IMU.StdevA value' )
    sdv(3) = IMU.StdevY;
end

[f1,c1] = size(y_out_data1);
uy = rand(f1,c1);
ny = zeros(f1,c1);

if strcmp(Noisealphadot_Flag,'NO') && strcmp(Noisexddot_Flag,'NO') && strcmp(Noiseyddot_Flag,'NO')
    out_data = y_out_data1;
    
elseif strcmp(Noisealphadot_Flag,'NO') && strcmp(Noisexddot_Flag,'NO') && strcmp(Noiseyddot_Flag,'YES')   
    ny(:,3) = icdf('Normal',uy(:,3),mu,sdv(3));
    out_data = y_out_data1 + ny;
    
elseif strcmp(Noisealphadot_Flag,'NO') && strcmp(Noisexddot_Flag,'YES') && strcmp(Noiseyddot_Flag,'NO')
    ny(:,2) = icdf('Normal',uy(:,2),mu,sdv(2));
    out_data = y_out_data1 + ny;
    
elseif strcmp(Noisealphadot_Flag,'NO') && strcmp(Noisexddot_Flag,'YES') && strcmp(Noiseyddot_Flag,'YES')
    ny(:,2) = icdf('Normal',uy(:,2),mu,sdv(2));
    ny(:,3) = icdf('Normal',uy(:,3),mu,sdv(3));
    out_data = y_out_data1 + ny;
    
elseif strcmp(Noisealphadot_Flag,'YES') && strcmp(Noisexddot_Flag,'NO') && strcmp(Noiseyddot_Flag,'NO')
    ny(:,1) = icdf('Normal',uy(:,1),mu,sdv(1));
    out_data = y_out_data1 + ny;

elseif strcmp(Noisealphadot_Flag,'YES') && strcmp(Noisexddot_Flag,'NO') && strcmp(Noiseyddot_Flag,'YES')
    ny(:,1) = icdf('Normal',uy(:,1),mu,sdv(1));
    ny(:,3) = icdf('Normal',uy(:,3),mu,sdv(3));
    out_data = y_out_data1 + ny;
    
elseif strcmp(Noisealphadot_Flag,'YES') && strcmp(Noisexddot_Flag,'YES') && strcmp(Noiseyddot_Flag,'NO')
    ny(:,1) = icdf('Normal',uy(:,1),mu,sdv(1));
    ny(:,2) = icdf('Normal',uy(:,2),mu,sdv(2));
    out_data = y_out_data1 + ny;
    
elseif strcmp(Noisealphadot_Flag,'YES') && strcmp(Noisexddot_Flag,'YES') && strcmp(Noiseyddot_Flag,'YES')
    ny(:,1) = icdf('Normal',uy(:,1),mu,sdv(1));
    ny(:,2) = icdf('Normal',uy(:,2),mu,sdv(2));
    ny(:,3) = icdf('Normal',uy(:,3),mu,sdv(3));
    out_data = y_out_data1 + ny;
    
else
    disp('WARNING: This Noisealphadot_Flag, Noisexddot_Flag and Noiseyddot_Flag combination does not exist, data will be sampled with no noise')
    out_data = y_out_data1;
end

end

