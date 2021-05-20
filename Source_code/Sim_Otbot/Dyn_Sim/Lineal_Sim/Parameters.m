%% System Matrices
Actc = [zeros(3,3),eye(3,3);
    zeros(3,3), zeros(3,3)];
Bctc = [zeros(3,3);eye(3,3)];

Cctc = diag(ones(6,1));
Dctc = zeros(6,3);

K1mat = diag([20,12.2,200]);
K2mat = diag([40,20.61,220]);
K = [K1mat,K2mat];

xs0=zeros(6,1);

%% Alternative control with LQR (in the nonlinear sistem we have oscil·lations)

% We have oscilations beacuse the poles were placed as complex (see the
% document eigvalsABK inside the folder tuning LQR)

%% Alternativa K when tuning with PP and using diferent poles
K1mat = diag([10,10,10]);
K2mat = diag([10,10,10]);
K = [K1mat,K2mat];
