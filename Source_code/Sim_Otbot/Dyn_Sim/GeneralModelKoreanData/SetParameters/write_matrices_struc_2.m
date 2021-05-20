load("m_struc")
load("dynamic_model_workspace")
load("forward_inverse_jacobians.mat")
%% Creating the missing symbolic variables

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

% Now MFIK matrix
MFIKmatrix = subs(MFIK_of_q,[l_1,l_2,r],[m.l_1,m.l_2,m.r]);

% Save matrices in a structure array
sm.MFIKmatrix = matlabFunction(MFIKmatrix);


%% Add the experssion of the Kinetics Energies

% TOTAL KYNETIC ENERGY OF THE SYSTEM
Texpr = subs(T,[I_b, I_p, I_a, I_t, l_1, l_2, m_b, m_w, m_p, x_G, y_G, x_F, y_F, r], [m.I_b, m.I_p, m.I_a, m.I_t, m.l_1, m.l_2, m.m_b, m.m_w, m.m_p, m.x_G, m.y_G, m.x_F, m.y_F, m.r]);
% Save kinetic energy expression in structure array
sm.Texpr = matlabFunction(Texpr);


%%%%%%%%%%%%%%%% ROTATION KYNETIC ENERGIES %%%%%%%%%%%%%%%%

% Rotation Energy of the chassis body
Trotbase_expr = subs(T_rot_b,[I_b, I_p, I_a, I_t, l_1, l_2, m_b, m_w, m_p, x_G, y_G, x_F, y_F, r], [m.I_b, m.I_p, m.I_a, m.I_t, m.l_1, m.l_2, m.m_b, m.m_w, m.m_p, m.x_G, m.y_G, m.x_F, m.y_F, m.r]);
%Save it
sm.Trotbase = matlabFunction(Trotbase_expr);

% Rotation Energy of the left wheel
Trot_leftwheel_expr = subs(T_rot_l, [I_b, I_p, I_a, I_t, l_1, l_2, m_b, m_w, m_p, x_G, y_G, x_F, y_F, r], [m.I_b, m.I_p, m.I_a, m.I_t, m.l_1, m.l_2, m.m_b, m.m_w, m.m_p, m.x_G, m.y_G, m.x_F, m.y_F, m.r]);
% Save it
sm.Trot_leftwheel_expr = matlabFunction(Trot_leftwheel_expr);

% Rotation Energy of the right wheel
Trot_rightwheel_expr = subs(T_rot_r, [I_b, I_p, I_a, I_t, l_1, l_2, m_b, m_w, m_p, x_G, y_G, x_F, y_F, r], [m.I_b, m.I_p, m.I_a, m.I_t, m.l_1, m.l_2, m.m_b, m.m_w, m.m_p, m.x_G, m.y_G, m.x_F, m.y_F, m.r]);
% Save it
sm.Trot_rightwheel_expr = matlabFunction(Trot_rightwheel_expr);

% Rotation Energy of the platform
Trot_platform_expr = subs(T_rot_p, [I_b, I_p, I_a, I_t, l_1, l_2, m_b, m_w, m_p, x_G, y_G, x_F, y_F, r], [m.I_b, m.I_p, m.I_a, m.I_t, m.l_1, m.l_2, m.m_b, m.m_w, m.m_p, m.x_G, m.y_G, m.x_F, m.y_F, m.r]);
% Save it
sm.Trot_platform_expr = matlabFunction(Trot_platform_expr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% TRANSLATION KYNETIC ENERGIES %%%%%%%%%%%%%%

% Translation Energy of the chassis body
Ttrabase_expr = subs(T_tra_b,[I_b, I_p, I_a, I_t, l_1, l_2, m_b, m_w, m_p, x_G, y_G, x_F, y_F, r], [m.I_b, m.I_p, m.I_a, m.I_t, m.l_1, m.l_2, m.m_b, m.m_w, m.m_p, m.x_G, m.y_G, m.x_F, m.y_F, m.r]);
%Save it
sm.Ttrabase = matlabFunction(Ttrabase_expr);

% Translation Energy of the left wheel
Ttra_leftwheel_expr = subs(T_tra_l,[I_b, I_p, I_a, I_t, l_1, l_2, m_b, m_w, m_p, x_G, y_G, x_F, y_F, r], [m.I_b, m.I_p, m.I_a, m.I_t, m.l_1, m.l_2, m.m_b, m.m_w, m.m_p, m.x_G, m.y_G, m.x_F, m.y_F, m.r]);
%Save it
sm.Ttra_leftwheel_expr = matlabFunction(Ttra_leftwheel_expr);

% Translation Energy of the right wheel
Ttra_rightwheel_expr = subs(T_tra_r,[I_b, I_p, I_a, I_t, l_1, l_2, m_b, m_w, m_p, x_G, y_G, x_F, y_F, r], [m.I_b, m.I_p, m.I_a, m.I_t, m.l_1, m.l_2, m.m_b, m.m_w, m.m_p, m.x_G, m.y_G, m.x_F, m.y_F, m.r]);
%Save it
sm.Ttra_rightwheel_expr = matlabFunction(Ttra_rightwheel_expr);

% Translation Energy of the platform
Ttra_platform_expr = subs(T_tra_p,[I_b, I_p, I_a, I_t, l_1, l_2, m_b, m_w, m_p, x_G, y_G, x_F, y_F, r], [m.I_b, m.I_p, m.I_a, m.I_t, m.l_1, m.l_2, m.m_b, m.m_w, m.m_p, m.x_G, m.y_G, m.x_F, m.y_F, m.r]);
%Save it
sm.Ttra_platform_expr = matlabFunction(Ttra_platform_expr);


%% Saveing only the structure object and deleting everything else
save('C:\Users\pereg\Documents\UPC\Master\TFMPractiques\MatlabWD\Sim_Otbot\Dyn_Sim\GeneralModelKoreanData\sm_struc.mat',"sm","sm")
clearvars
close all
clc
