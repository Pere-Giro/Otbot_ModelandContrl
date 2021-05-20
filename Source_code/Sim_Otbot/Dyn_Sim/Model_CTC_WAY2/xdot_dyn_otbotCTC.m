function [ xdot ] = xdot_dyn_otbotCTC(t, xs, sm, u, u_Flag, cp, goal_Flag)
%XDOT_HEXAPOLE Summary of this function goes here
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
% varphi_dot_r = xs(10);
% varphi_dot_l = xs(11);
varphi_dot_p = xs(12);

MIIKs = sm.MIIKmatrix(alpha,varphi_p);
MIIKdots = sm.MIIKdotamatrix(alpha,alpha_dot,varphi_p,varphi_dot_p);




qdot = xs(7:12);

switch u_Flag
    case 'CTE'
        qdotdot(1:3,1) = u;
        qdotdot(4:6,1) = MIIKdots*qdot(1:3,1) + MIIKs*u;
    case 'VAR'
        u_f= u_function(t,u);
        qdotdot(1:3,1) = u_f;
        qdotdot(4:6,1) = MIIKdots*qdot(1:3,1) + MIIKs*u_f;
        
    case 'CTC_LQR'
        switch goal_Flag
            case 'FIX'
                pdiff(1:2,1) = cp.pss(1:2,1)-xs(1:2,1);
                %pdiff(3,1) = angdiff(xs(3,1), cp.pss(3,1));
                pdiff(3,1) = cp.pss(3,1)-xs(3,1);
                pdiff(4:6,1) = cp.pss(4:6,1)-xs(7:9,1);

                u = cp.K*pdiff;       
                qdotdot(1:3,1) = u;
                qdotdot(4:6,1) = MIIKdots*qdot(1:3,1) + MIIKs*u;
                
            case 'CIRCLE'
                [xss,yss,alphass,xdotss,ydotss,alphadotss,xdotdotss,ydotdotss,alphadotdotss] = TD_circle_of_time_XYA(t);
                cp.pss = [xss;yss;alphass;xdotss;ydotss;alphadotss];
                pd_dotdot = [xdotdotss; ydotdotss; alphadotdotss];
                
                pdiff(1:2,1) = cp.pss(1:2,1)-xs(1:2,1);
                %pdiff(3,1) = angdiff(xs(3,1), cp.pss(3,1));
                pdiff(3,1) = cp.pss(3,1)-xs(3,1);
                pdiff(4:6,1) = cp.pss(4:6,1)-xs(7:9,1);

                u = pd_dotdot + cp.K*pdiff;
                qdotdot(1:3,1) = u;
                qdotdot(4:6,1) = MIIKdots*qdot(1:3,1) + MIIKs*u;
                
            case 'POLILINE'
                [xss,yss,alphass,xdotss,ydotss,alphadotss,xdotdotss,ydotdotss,alphadotdotss] = TD_polyline_of_time(t);
                cp.pss = [xss;yss;alphass;xdotss;ydotss;alphadotss];
                pd_dotdot = [xdotdotss; ydotdotss; alphadotdotss];
                
                pdiff(1:2,1) = cp.pss(1:2,1)-xs(1:2,1);
                %pdiff(3,1) = angdiff(xs(3,1), cp.pss(3,1));
                pdiff(3,1) = cp.pss(3,1)-xs(3,1);
                pdiff(4:6,1) = cp.pss(4:6,1)-xs(7:9,1);

                u = pd_dotdot + cp.K*pdiff;
                qdotdot(1:3,1) = u;
                qdotdot(4:6,1) = MIIKdots*qdot(1:3,1) + MIIKs*u;
                
            otherwise
                disp('This goal_Flag does not exist (xdot_dyn_otbotCTC)')
                disp('Simulating with a FIX point X=2 Y=2 Alpha = 0')
                cp.pss = [2,2,0,0,0,0]';
                
                pdiff(1:2,1) = cp.pss(1:2,1)-xs(1:2,1);
                %pdiff(3,1) = angdiff(xs(3,1), cp.pss(3,1));
                pdiff(3,1) = cp.pss(3,1)-xs(3,1);
                pdiff(4:6,1) = cp.pss(4:6,1)-xs(7:9,1);

                u = cp.K*pdiff;       
                qdotdot(1:3,1) = u;
                qdotdot(4:6,1) = MIIKdots*qdot(1:3,1) + MIIKs*u;
        end
        
    case 'CTC_PP' %Note that this section is a copy exacly of the CTC_LQR
        switch goal_Flag
            case 'FIX'
                pdiff(1:2,1) = cp.pss(1:2,1)-xs(1:2,1);
                %pdiff(3,1) = angdiff(xs(3,1), cp.pss(3,1));
                pdiff(3,1) = cp.pss(3,1)-xs(3,1);
                pdiff(4:6,1) = cp.pss(4:6,1)-xs(7:9,1);

                u = cp.K*pdiff;
                qdotdot(1:3,1) = u;
                qdotdot(4:6,1) = MIIKdots*qdot(1:3,1) + MIIKs*u;
                
            case 'CIRCLE'
                [xss,yss,alphass,xdotss,ydotss,alphadotss,xdotdotss,ydotdotss,alphadotdotss] = TD_circle_of_time_XYA(t);
                cp.pss = [xss;yss;alphass;xdotss;ydotss;alphadotss];
                pd_dotdot = [xdotdotss; ydotdotss; alphadotdotss];
                
                pdiff(1:2,1) = cp.pss(1:2,1)-xs(1:2,1);
                %pdiff(3,1) = angdiff(xs(3,1), cp.pss(3,1));
                pdiff(3,1) = cp.pss(3,1)-xs(3,1);
                pdiff(4:6,1) = cp.pss(4:6,1)-xs(7:9,1);

                u = pd_dotdot + cp.K*pdiff;
                qdotdot(1:3,1) = u;
                qdotdot(4:6,1) = MIIKdots*qdot(1:3,1) + MIIKs*u;
                
            case 'POLILINE'
                [xss,yss,alphass,xdotss,ydotss,alphadotss,xdotdotss,ydotdotss,alphadotdotss] = TD_polyline_of_time(t);
                cp.pss = [xss;yss;alphass;xdotss;ydotss;alphadotss];
                pd_dotdot = [xdotdotss; ydotdotss; alphadotdotss];
                
                pdiff(1:2,1) = cp.pss(1:2,1)-xs(1:2,1);
                %pdiff(3,1) = angdiff(xs(3,1), cp.pss(3,1));
                pdiff(3,1) = cp.pss(3,1)-xs(3,1);
                pdiff(4:6,1) = cp.pss(4:6,1)-xs(7:9,1);

                u = pd_dotdot + cp.K*pdiff;
                qdotdot(1:3,1) = u;
                qdotdot(4:6,1) = MIIKdots*qdot(1:3,1) + MIIKs*u;
                
            otherwise
                disp('This goal_Flag does not exist (xdot_dyn_otbotCTC)')
                disp('Simulating with a FIX point X=2 Y=2 Alpha = 0')
                cp.pss = [2,2,0,0,0,0]';
                
                pdiff(1:2,1) = cp.pss(1:2,1)-xs(1:2,1);
                %pdiff(3,1) = angdiff(xs(3,1), cp.pss(3,1));
                pdiff(3,1) = cp.pss(3,1)-xs(3,1);
                pdiff(4:6,1) = cp.pss(4:6,1)-xs(7:9,1);

                u = cp.K*pdiff;       
                qdotdot(1:3,1) = u;
                qdotdot(4:6,1) = MIIKdots*qdot(1:3,1) + MIIKs*u;
        end
    otherwise
        disp('WARNING THIS FLAG DOES NOT EXIST (xdot_dyn_otbotCTC)')
        disp('Simulation with CTE torques will be launched')
        qdotdot(1:3,1) = u;
        qdotdot(4:6,1) = MIIKdots*qdot(1:3,1) + MIIKs*u;
end
xdot=[qdot; qdotdot];
end





