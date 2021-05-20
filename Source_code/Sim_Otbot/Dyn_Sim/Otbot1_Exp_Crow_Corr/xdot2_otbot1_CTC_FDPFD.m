function [ xdot ] = xdot2_otbot1_CTC_FDPFD(t, xs, m, sm, u, u_Flag, cp, goal_Flag, FrictFlag, DistFlag, u_IntFlag, Int_Vec,uactionts)
% XDOT_HEXAPOLE Summary of this function goes here
%   differential equation (non-linearized), that modelises our hexapole
%   system

% x = xs(1);
% y = xs(2);
alpha = xs(3);
% varphi_r = xs(4);
% varphi_l = xs(5);
varphi_p = xs(6);

% x_dot = xs(7);
% y_dot = xs(8);
alpha_dot = xs(9);
varphi_dot_r = xs(10);
varphi_dot_l = xs(11);
varphi_dot_p = xs(12);

% Setting fricctions
switch FrictFlag
    case 'YES'
        tau_friction(1:3,1) = zeros(3,1);
        tau_friction(4,1) = -m.b_frict(1,1)*varphi_dot_r; % tau_friction_right
        tau_friction(5,1) = -m.b_frict(2,1)*varphi_dot_l; % tau_friction_left
        tau_friction(6,1) = -m.b_frict(3,1)*varphi_dot_p; % tau_friction_pivot
    case 'NO'
        tau_friction = zeros(6,1);
    otherwise
        disp('Warning: This FrictFlag does not exist. Sistem will be simulated without friction')
        tau_friction = zeros(6,1);
end

% Setting Distrurbances
switch DistFlag
    case 'YES'
        D_vec = Dist_funct(t);
    case 'NO'
        D_vec = zeros(6,1);
    otherwise
        disp('Warning: This DistFlag does not exist. Sistim will be simulated without disturbances')
        D_vec = zeros(6,1);
end
     
MIIKs = sm.MIIKmatrix(alpha,varphi_p);
M_bars2 = sm.M_barmatrix2(alpha,varphi_p);
C_bars2 = sm.C_barmatrix2(alpha,alpha_dot,varphi_p,varphi_dot_p);
MIIKdots = sm.MIIKdotmatrix(alpha,alpha_dot,varphi_p,varphi_dot_p);
Dmats = sm.Deltamatrix(alpha,varphi_p);

Msys2 = [M_bars2, zeros(3);
         -MIIKs, eye(3)];

qdot = xs(7:12);

switch u_Flag
    case 'CTE'
        qdotdot = pinv(Msys2)*[u + Dmats.'*tau_friction + Dmats.'*D_vec - C_bars2*qdot(1:3,1); MIIKdots*qdot(1:3,1)]; 
        
    case 'VAR'
        u_f = u_function(t,uactionts);
        qdotdot = pinv(Msys2)*[u_f + Dmats.'*tau_friction + Dmats.'*D_vec - C_bars2*qdot(1:3,1); MIIKdots*qdot(1:3,1)];
                
    case 'CTC_LQR'
        switch goal_Flag
            case 'FIX'
                pdiff(1:2,1) = cp.pss(1:2,1)-xs(1:2,1);
                %pdiff(3,1) = angdiff(xs(3,1), cp.pss(3,1));
                pdiff(3,1) = cp.pss(3,1)-xs(3,1);
                pdiff(4:6,1) = cp.pss(4:6,1)-xs(7:9,1);
                
                u_2 = cp.K*pdiff;
                u = M_bars2*u_2 + C_bars2*qdot(1:3,1);        
                qdotdot = pinv(Msys2)*[u + Dmats.'*tau_friction + Dmats.'*D_vec - C_bars2*qdot(1:3,1); MIIKdots*qdot(1:3,1)];                
                
            case 'CIRCLE'
                [xss,yss,alphass,xdotss,ydotss,alphadotss,xdotdotss,ydotdotss,alphadotdotss] = TD_circle_of_time_XYA(t);
                cp.pss = [xss;yss;alphass;xdotss;ydotss;alphadotss];
                pd_dotdot = [xdotdotss; ydotdotss; alphadotdotss];
                
                pdiff(1:2,1) = cp.pss(1:2,1)-xs(1:2,1);
                %pdiff(3,1) = angdiff(xs(3,1), cp.pss(3,1));
                pdiff(3,1) = cp.pss(3,1)-xs(3,1);
                pdiff(4:6,1) = cp.pss(4:6,1)-xs(7:9,1);

                u_2 = pd_dotdot + cp.K*pdiff;
                u = M_bars2*u_2 + C_bars2*qdot(1:3,1);        
                qdotdot = pinv(Msys2)*[u + Dmats.'*tau_friction + Dmats.'*D_vec - C_bars2*qdot(1:3,1); MIIKdots*qdot(1:3,1)];
                                
            case 'POLILINE'
                [xss,yss,alphass,xdotss,ydotss,alphadotss,xdotdotss,ydotdotss,alphadotdotss] = TD_polyline_of_time(t);
                cp.pss = [xss;yss;alphass;xdotss;ydotss;alphadotss];
                pd_dotdot = [xdotdotss; ydotdotss; alphadotdotss];
                
                pdiff(1:2,1) = cp.pss(1:2,1)-xs(1:2,1);
                %pdiff(3,1) = angdiff(xs(3,1), cp.pss(3,1));
                pdiff(3,1) = cp.pss(3,1)-xs(3,1);
                pdiff(4:6,1) = cp.pss(4:6,1)-xs(7:9,1);

                u_2 = pd_dotdot + cp.K*pdiff;
                u = M_bars2*u_2 + C_bars2*qdot(1:3,1);        
                qdotdot = pinv(Msys2)*[u + Dmats.'*tau_friction + Dmats.'*D_vec - C_bars2*qdot(1:3,1); MIIKdots*qdot(1:3,1)];
                                
            otherwise
                disp('This goal_Flag does not exist (xdot_dyn_otbotCTC)')
                disp('Simulating with a FIX point X=2 Y=2 Alpha = 0')
                cp.pss = [2,2,0,0,0,0]';
                
                pdiff(1:2,1) = cp.pss(1:2,1)-xs(1:2,1);
                %pdiff(3,1) = angdiff(xs(3,1), cp.pss(3,1));
                pdiff(3,1) = cp.pss(3,1)-xs(3,1);
                pdiff(4:6,1) = cp.pss(4:6,1)-xs(7:9,1);

                u_2 = cp.K*pdiff;
                u = M_bars2*u_2 + C_bars2*qdot(1:3,1);        
                qdotdot = pinv(Msys2)*[u + Dmats.'*tau_friction + Dmats.'*D_vec - C_bars2*qdot(1:3,1); MIIKdots*qdot(1:3,1)];
                
        end
        
    case 'CTC_PP' % This secction used to be a copy of LQR but now it is included the integrative term
        switch goal_Flag
            case 'FIX'
                pdiff(1:2,1) = cp.pss(1:2,1)-xs(1:2,1);
                %pdiff(3,1) = angdiff(xs(3,1), cp.pss(3,1));
                pdiff(3,1) = cp.pss(3,1)-xs(3,1);
                pdiff(4:6,1) = cp.pss(4:6,1)-xs(7:9,1);
                
                if strcmp(u_IntFlag,'YES')
                    I_k = Int_Vec(1,1:3)' + (t-Int_Vec(1,7))/2*(Int_Vec(1,4:6)' + pdiff(1:3,1));
                    u_2 = cp.K*[I_k;pdiff];
                    u = M_bars2*u_2 + C_bars2*qdot(1:3,1);        
                    qdotdot = pinv(Msys2)*[u + Dmats.'*tau_friction + Dmats.'*D_vec - C_bars2*qdot(1:3,1); MIIKdots*qdot(1:3,1)];
                    
                    Int_Vec(1,1:3) = I_k';          % Updating Integral value
                    Int_Vec(1,4:6) = pdiff(1:3,1)'; % Updating error value
                    Int_Vec(1,7) = t;               % Updating time value
                elseif strcmp(u_IntFlag,'NO')
                    u_2 = cp.K*pdiff;
                    u = M_bars2*u_2 + C_bars2*qdot(1:3,1);        
                    qdotdot = pinv(Msys2)*[u + Dmats.'*tau_friction + Dmats.'*D_vec - C_bars2*qdot(1:3,1); MIIKdots*qdot(1:3,1)];
                else
                    disp('This u_IntFlag does not exist runing the simulation with standart PD controller')
                    u_2 = cp.K*pdiff;
                    u = M_bars2*u_2 + C_bars2*qdot(1:3,1);        
                    qdotdot = pinv(Msys2)*[u + Dmats.'*tau_friction + Dmats.'*D_vec - C_bars2*qdot(1:3,1); MIIKdots*qdot(1:3,1)];
                end
                
%                 u_2 = cp.K*pdiff;
%                 u = (MFIKs.')*(M_bars*u_2 + C_bars*qdot(1:3,1));        
%                 qdotdot = [eye(6),zeros(6,3)]*pinv([Ms,Js.'; Js, zeros(3,3)])*[Es*u + tau_friction + D_vec - Cs*qdot; -Jdots*qdot];
                                
            case 'CIRCLE'
                [xss,yss,alphass,xdotss,ydotss,alphadotss,xdotdotss,ydotdotss,alphadotdotss] = TD_circle_of_time_XYA(t);
                cp.pss = [xss;yss;alphass;xdotss;ydotss;alphadotss];
                pd_dotdot = [xdotdotss; ydotdotss; alphadotdotss];
                
                pdiff(1:2,1) = cp.pss(1:2,1)-xs(1:2,1);
                %pdiff(3,1) = angdiff(xs(3,1), cp.pss(3,1));
                pdiff(3,1) = cp.pss(3,1)-xs(3,1);
                pdiff(4:6,1) = cp.pss(4:6,1)-xs(7:9,1);
                
                if strcmp(u_IntFlag,'YES')
                    I_k = Int_Vec(1,1:3)' + (t-Int_Vec(1,7))/2*(Int_Vec(1,4:6)' + pdiff(1:3,1));
                    u_2 = pd_dotdot + cp.K*[I_k;pdiff];
                    u = M_bars2*u_2 + C_bars2*qdot(1:3,1);        
                    qdotdot = pinv(Msys2)*[u + Dmats.'*tau_friction + Dmats.'*D_vec - C_bars2*qdot(1:3,1); MIIKdots*qdot(1:3,1)];
                    
                    Int_Vec(1,1:3) = I_k';          % Updating Integral value
                    Int_Vec(1,4:6) = pdiff(1:3,1)'; % Updating error value
                    Int_Vec(1,7) = t;               % Updating time value
                elseif strcmp(u_IntFlag,'NO')
                    u_2 = pd_dotdot + cp.K*pdiff;
                    u = M_bars2*u_2 + C_bars2*qdot(1:3,1);        
                    qdotdot = pinv(Msys2)*[u + Dmats.'*tau_friction + Dmats.'*D_vec - C_bars2*qdot(1:3,1); MIIKdots*qdot(1:3,1)];
                else
                    disp('This u_IntFlag does not exist runing the simulation with standart PD controller')
                    u_2 = pd_dotdot + cp.K*pdiff;
                    u = M_bars2*u_2 + C_bars2*qdot(1:3,1);        
                    qdotdot = pinv(Msys2)*[u + Dmats.'*tau_friction + Dmats.'*D_vec - C_bars2*qdot(1:3,1); MIIKdots*qdot(1:3,1)];
                end
                
%                 u_2 = pd_dotdot + cp.K*pdiff;
%                 u = (MFIKs.')*(M_bars*u_2 + C_bars*qdot(1:3,1));        
%                 qdotdot = [eye(6),zeros(6,3)]*pinv([Ms,Js.'; Js, zeros(3,3)])*[Es*u + tau_friction + D_vec - Cs*qdot; -Jdots*qdot];
                                
            case 'POLILINE'
                [xss,yss,alphass,xdotss,ydotss,alphadotss,xdotdotss,ydotdotss,alphadotdotss] = TD_polyline_of_time(t);
                cp.pss = [xss;yss;alphass;xdotss;ydotss;alphadotss];
                pd_dotdot = [xdotdotss; ydotdotss; alphadotdotss];
                
                pdiff(1:2,1) = cp.pss(1:2,1)-xs(1:2,1);
                %pdiff(3,1) = angdiff(xs(3,1), cp.pss(3,1));
                pdiff(3,1) = cp.pss(3,1)-xs(3,1);
                pdiff(4:6,1) = cp.pss(4:6,1)-xs(7:9,1);
                
                if strcmp(u_IntFlag,'YES')
                    I_k = Int_Vec(1,1:3)' + (t-Int_Vec(1,7))/2*(Int_Vec(1,4:6)' + pdiff(1:3,1));
                    u_2 = pd_dotdot + cp.K*[I_k;pdiff];
                    u = M_bars2*u_2 + C_bars2*qdot(1:3,1);        
                    qdotdot = pinv(Msys2)*[u + Dmats.'*tau_friction + Dmats.'*D_vec - C_bars2*qdot(1:3,1); MIIKdots*qdot(1:3,1)];
                    
                    Int_Vec(1,1:3) = I_k';          % Updating Integral value
                    Int_Vec(1,4:6) = pdiff(1:3,1)'; % Updating error value
                    Int_Vec(1,7) = t;               % Updating time value
                elseif strcmp(u_IntFlag,'NO')
                    u_2 = pd_dotdot + cp.K*pdiff;
                    u = M_bars2*u_2 + C_bars2*qdot(1:3,1);        
                    qdotdot = pinv(Msys2)*[u + Dmats.'*tau_friction + Dmats.'*D_vec - C_bars2*qdot(1:3,1); MIIKdots*qdot(1:3,1)];
                else
                    disp('This u_IntFlag does not exist runing the simulation with standart PD controller')
                    u_2 = pd_dotdot + cp.K*pdiff;
                    u = M_bars2*u_2 + C_bars2*qdot(1:3,1);        
                    qdotdot = pinv(Msys2)*[u + Dmats.'*tau_friction + Dmats.'*D_vec - C_bars2*qdot(1:3,1); MIIKdots*qdot(1:3,1)];
                end

%                 u_2 = pd_dotdot + cp.K*pdiff;
%                 u = (MFIKs.')*(M_bars*u_2 + C_bars*qdot(1:3,1));        
%                 qdotdot = [eye(6),zeros(6,3)]*pinv([Ms,Js.'; Js, zeros(3,3)])*[Es*u + tau_friction + D_vec - Cs*qdot; -Jdots*qdot];
                                
            otherwise
                disp('This goal_Flag does not exist (xdot_dyn_otbotCTC)')
                disp('Simulating with a FIX point X=2 Y=2 Alpha = 0 with standart PD CTC')
                cp.pss = [2,2,0,0,0,0]';
                
                pdiff(1:2,1) = cp.pss(1:2,1)-xs(1:2,1);
                %pdiff(3,1) = angdiff(xs(3,1), cp.pss(3,1));
                pdiff(3,1) = cp.pss(3,1)-xs(3,1);
                pdiff(4:6,1) = cp.pss(4:6,1)-xs(7:9,1);

                u_2 = cp.K*pdiff;
                u = M_bars2*u_2 + C_bars2*qdot(1:3,1);        
                qdotdot = pinv(Msys2)*[u + Dmats.'*tau_friction + Dmats.'*D_vec - C_bars2*qdot(1:3,1); MIIKdots*qdot(1:3,1)];
                                
        end
        
    otherwise
        disp('WARNING THIS FLAG DOES NOT EXIST (xdot_dyn_otbotCTC)')
        disp('Simulation with CTE torques will be launched')
        qdotdot = pinv(Msys2)*[u + Dmats.'*tau_friction + Dmats.'*D_vec - C_bars2*qdot(1:3,1); MIIKdots*qdot(1:3,1)];
        
end
xdot=[qdot; qdotdot];
end





