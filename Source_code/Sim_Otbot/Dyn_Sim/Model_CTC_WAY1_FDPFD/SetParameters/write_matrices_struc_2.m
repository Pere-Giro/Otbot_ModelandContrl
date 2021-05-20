load("m_struc")
load("dynamic_model_workspace")
load("M_bar.mat")
load("C_bar.mat")
load("forward_inverse_jacobians.mat")
load("Lambmat.mat")

%% Create missing simbolic variables

syms r

%% Seting up

% Now we will be seting up this matrices in a apropiate way in order to use
% them for the simulations.

% First J matrix
Jmatrix = subs(J_of_q,[l_1,l_2,r],[m.l_1,m.l_2,m.r]);

% Save matrices in a structure array
sm.Jmatrix = matlabFunction(Jmatrix);

% Now Jdot matrix
Jdotmatrix = subs(Jdot,[l_1,l_2,r],[m.l_1,m.l_2,m.r]);

% Save matrices in a structure array
sm.Jdotmatrix = matlabFunction(Jdotmatrix);

% Now Mass matrix
Mmatrix = subs(M,[I_b, I_p, I_a, I_t, l_1, l_2, m_b, m_w, m_p, x_G, y_G, x_F, y_F, r], [m.I_b, m.I_p, m.I_a, m.I_t, m.l_1, m.l_2, m.m_b, m.m_w, m.m_p, m.x_G, m.y_G, m.x_F, m.y_F, m.r]);

% Save matrices in a structure array
sm.Mmatrix = matlabFunction(Mmatrix);

%Now we include E matrix also
sm.Ematrix = Ematrix;

% Now Cmat
Cmatrix = subs(Cmat,[I_b, I_p, I_a, I_t, l_1, l_2, m_b, m_w, m_p, x_G, y_G, x_F, y_F, r], [m.I_b, m.I_p, m.I_a, m.I_t, m.l_1, m.l_2, m.m_b, m.m_w, m.m_p, m.x_G, m.y_G, m.x_F, m.y_F, m.r]);

% Save matrices in a structure array
sm.Cmatrix = matlabFunction(Cmatrix);

%% For this taskSpace model we need to use task space matrices

% Let's set-up M_bar and C_bar and even MIIK

% First of all C_bar
C_barmatrix = subs(C_bar,[I_b, I_p, I_a, I_t, l_1, l_2, m_b, m_w, m_p, x_G, y_G, x_F, y_F, r], [m.I_b, m.I_p, m.I_a, m.I_t, m.l_1, m.l_2, m.m_b, m.m_w, m.m_p, m.x_G, m.y_G, m.x_F, m.y_F, m.r]);

% Save matrices in a structure array
sm.C_barmatrix = matlabFunction(C_barmatrix);

% Then M_bar
M_barmatrix = subs(M_bar,[I_b, I_p, I_a, I_t, l_1, l_2, m_b, m_w, m_p, x_G, y_G, x_F, y_F, r], [m.I_b, m.I_p, m.I_a, m.I_t, m.l_1, m.l_2, m.m_b, m.m_w, m.m_p, m.x_G, m.y_G, m.x_F, m.y_F, m.r]);

% Save matrices in a structure array
sm.M_barmatrix = matlabFunction(M_barmatrix);

% MIIK matrix
MIIKmatrix = subs(MIIK_of_q,[l_1,l_2,r],[m.l_1,m.l_2,m.r]);

% Save matrices in a structure array
sm.MIIKmatrix = matlabFunction(MIIKmatrix);

% MFIK matrix
MFIKmatrix = subs(MFIK_of_q,[l_1,l_2,r],[m.l_1,m.l_2,m.r]);

% Save matrices in a structure array
sm.MFIKmatrix = matlabFunction(MFIKmatrix);

% Lambda Matrix
Lambdamatrix = subs(Lambda,[l_1,l_2,r],[m.l_1,m.l_2,m.r]);

% Save matrix in a structure array
sm.Lambdamatrix = matlabFunction(Lambdamatrix);

%% Add the experssion of the Kinetic Energy
Texpr = subs(T,[I_b, I_p, I_a, I_t, l_1, l_2, m_b, m_w, m_p, x_G, y_G, x_F, y_F, r], [m.I_b, m.I_p, m.I_a, m.I_t, m.l_1, m.l_2, m.m_b, m.m_w, m.m_p, m.x_G, m.y_G, m.x_F, m.y_F, m.r]);
% Save kinetic energy expression in structure array
sm.Texpr = matlabFunction(Texpr);

%% Saveing only the structure object and deleting everything else
save('C:\Users\pereg\Documents\UPC\Master\TFMPractiques\MatlabWD\Sim_Otbot\Dyn_Sim\Model_CTC_WAY1_DF\sm_struc.mat',"sm","sm")
clearvars
close all
clc
