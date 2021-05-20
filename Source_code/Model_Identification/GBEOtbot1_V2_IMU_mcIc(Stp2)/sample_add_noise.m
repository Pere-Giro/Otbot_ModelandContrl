%% This script adds noise to the sample data

load("Otbot1_Data.mat")

% Normal noise parameters for the control inputs:
mu = 0; % mean

[f1,c1] = size(y_out_data);
[f2,c2] = size(u_out_data);

uy = rand(f1,c1);
uu = rand(f2,c2);

ny = zeros(f1,c1);
nu = zeros(f2,c2);

for j=1:c1
    max1 = max(y_out_data(:,j));
    min1 = min(y_out_data(:,j));
    if abs(max1)>1e-3 || abs(min1)>1e-3
        sdv1 = abs(max1-min1)/40;
        ny(:,j) = icdf('Normal',uy(:,j),mu,sdv1);
    end
end

for j=1:c2
    max2 = max(u_out_data(:,j));
    min2 = min(u_out_data(:,j));
    if abs(max2)>1e-3 || abs(min2)>1e-3
        sdv2 = abs(max1-min1)/40;
        if sdv2 < 1e-3
            sdv2 = 0.1;
        end
        nu(:,j) = icdf('Normal',uu(:,j),mu,sdv2);
    end
end

y_out_data = y_out_data + ny;
u_out_data = u_out_data + nu;

save("Otbot1_DataNoise.mat")
