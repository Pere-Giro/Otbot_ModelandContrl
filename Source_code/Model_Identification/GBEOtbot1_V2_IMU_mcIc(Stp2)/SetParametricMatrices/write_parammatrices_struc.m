load("M_bar_V2.mat")
load("C_bar_V2.mat")
load("forward_inverse_jacobians.mat")
load("MIIKdotmat.mat")
load("Lambmat.mat")
load("Deltamat.mat")
%% Creating the missing symbolic variables

syms r

%% Seting up

% Now we will be seting up this matrices in a apropiate way in order to use
% them for the simulations.

% Save matrices in a structure array
sm.C_barmatrix2 = matlabFunction(C_bar_V2);

% Save matrices in a structure array
sm.M_barmatrix2 = matlabFunction(M_bar_V2);

% Save matrices in a structure array
sm.MIIKmatrix = matlabFunction(MIIK_of_q);

% Save matrices in a structure array
sm.MIIKdotmatrix = matlabFunction(MIIKdot);

% Save matrices in a structure array
sm.MFIKmatrix = matlabFunction(MFIK_of_q);

% Save matrix in a structure array
sm.Lambdamatrix = matlabFunction(Lambda);

% Save matrix in a structure array
sm.Deltamatrix = matlabFunction(Delta);

%% Saveing only the structure object and deleting everything else
% save('C:\Users\pereg\Documents\UPC\Master\TFMPractiques\MatlabWD\Model_Identification\GBEOtbot1_V2_IMU_mcIc\sm_struc.mat',"sm","sm")
Upath = userpath;
savedirsp1 = strcat(Upath,'\Model_Identification\GBEOtbot1_V2_IMU_mcIc(Stp2)\');
save(strcat(savedirsp1,'sm_struc.mat'),"sm")
clearvars
close all
clc
