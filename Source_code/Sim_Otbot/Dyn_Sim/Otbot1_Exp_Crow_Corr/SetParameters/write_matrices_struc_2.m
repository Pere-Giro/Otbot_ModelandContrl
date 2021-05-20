load("m_struc")
load("dynamic_model_workspace")
load("M_bar_V2.mat")
load("C_bar_V2.mat")
load("forward_inverse_jacobians.mat")
load("Lambmat.mat")
load("Deltamat.mat")
load("MIIKdotmat.mat")

%% Create missing simbolic variables

syms r

%% Seting up

% Now we will be seting up this matrices in a apropiate way in order to use
% them for the simulations.

% First of all C_bar2
C_barmatrix2 = subs(C_bar_V2,[I_b, I_p, I_a, I_t, l_1, l_2, m_b, m_w, m_p, x_G, y_G, x_F, y_F, r], [m.I_b, m.I_p, m.I_a, m.I_t, m.l_1, m.l_2, m.m_b, m.m_w, m.m_p, m.x_G, m.y_G, m.x_F, m.y_F, m.r]);

% Save matrices in a structure array
sm.C_barmatrix2 = matlabFunction(C_barmatrix2);

% Then M_bar2
M_barmatrix2 = subs(M_bar_V2,[I_b, I_p, I_a, I_t, l_1, l_2, m_b, m_w, m_p, x_G, y_G, x_F, y_F, r], [m.I_b, m.I_p, m.I_a, m.I_t, m.l_1, m.l_2, m.m_b, m.m_w, m.m_p, m.x_G, m.y_G, m.x_F, m.y_F, m.r]);

% Save matrices in a structure array
sm.M_barmatrix2 = matlabFunction(M_barmatrix2);

% MIIK matrix
MIIKmatrix = subs(MIIK_of_q,[l_1,l_2,r],[m.l_1,m.l_2,m.r]);

% Save matrices in a structure array
sm.MIIKmatrix = matlabFunction(MIIKmatrix);

% MIIKdot matrix
MIIKdotmatrix = subs(MIIKdot,[l_1,l_2,r],[m.l_1,m.l_2,m.r]);

% Save matrices in a structure array
sm.MIIKdotmatrix = matlabFunction(MIIKdotmatrix);

% MFIK matrix
MFIKmatrix = subs(MFIK_of_q,[l_1,l_2,r],[m.l_1,m.l_2,m.r]);

% Save matrices in a structure array
sm.MFIKmatrix = matlabFunction(MFIKmatrix);

% Lambda Matrix
Lambdamatrix = subs(Lambda,[l_1,l_2,r],[m.l_1,m.l_2,m.r]);

% Save matrix in a structure array
sm.Lambdamatrix = matlabFunction(Lambdamatrix);

% Lambda Matrix
Deltamatrix = subs(Delta,[l_1,l_2,r],[m.l_1,m.l_2,m.r]);

% Save matrix in a structure array
sm.Deltamatrix = matlabFunction(Deltamatrix);


%% Add the experssion of the Kinetic Energy
Texpr = subs(T,[I_b, I_p, I_a, I_t, l_1, l_2, m_b, m_w, m_p, x_G, y_G, x_F, y_F, r], [m.I_b, m.I_p, m.I_a, m.I_t, m.l_1, m.l_2, m.m_b, m.m_w, m.m_p, m.x_G, m.y_G, m.x_F, m.y_F, m.r]);
% Save kinetic energy expression in structure array
sm.Texpr = matlabFunction(Texpr);

%% Saveing only the structure object and deleting everything else
save('C:\Users\pereg\Documents\UPC\Master\TFMPractiques\MatlabWD\Sim_Otbot\Dyn_Sim\Otbot1_Crow_Corr\sm_struc.mat',"sm")
clearvars
close all
clc
