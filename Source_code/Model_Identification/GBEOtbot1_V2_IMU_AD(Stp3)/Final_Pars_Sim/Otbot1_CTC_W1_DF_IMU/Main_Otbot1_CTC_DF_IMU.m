clearvars
close all
clc
%% WARNING
% This simulation uses always a zoh to give the input action to the
% system, meaning that the actions are not updated at the same frequancy
% that the evolution of the system itself

%% Load the parameters in the model
addpath('SetParameters') % Adding set parameters to search path
otbot1_parameters_1
rmpath('SetParameters') % removing set parameters of the search path
%% Flags
load('Flags_Data.mat'); % Loading the file with the flags
               
FrictFlag = Flags.FrictFlag; % Flag to define if we want to add friction in our model or not
               % (FrictFlag == YES)> Our model will have freictions
               % (FrictFlag == NO)> Our model will not contain frictions
               
DistFlag = Flags.DistFlag; % Flag to define if we want to add non modelled disturbances during the simulation
               % (DistFlag == YES)> Non modelled disturbances will be added 
               % (DistFlag == NO)> Non modelled disturbances will NOT be added
               
TrackFlag = Flags.TrackFlag; % Flag to define if you want to track the path of each center of mass
               % (TrackFlag == NO)> plots without any paths of center of
               % mass plotted              
               % (TrackFlag == CB)> plots the path of the center of mass of
               % the chassis body               
               % (TrackFlag == PF)> plots the path of the center of mass of
               % the platform body              
               % (TrackFlag == BOTH)> plots the path of both centers of
               % mass
               
TargetFlag = Flags.TargetFlag; % Flag to choose if you want to see the desired goal in the video simulation
                    % (Only avaliable if u_Flag is CTC_LQR or CTC_PP)
                % (TargetFlag == YES)> A red target will be displayed in the
                % simulation
                % (TargetFlag == NO)> There won't be a target in the video
                % simulation
               
UvecFlag = Flags.UvecFlag; % Flag to choose if you want to compute and display the plots of action vector u
                % (UvecFlag == YES)> action vector u plots will be displayed
                % (UvecFlag == NO)> action vector u plots will not be displayed
                           
KinEnFlag = Flags.KinEnFlag; % Flag to choose if Kinetic energic plot must be displayed or not
               % (KinEnFlag == YES)> Outputs a plot of kinetic energy
               % (KinEnFlag == NO )> No output plot of kinetic energy
               
Jdot_qdotFlag = Flags.Jdot_qdotFlag; % Flag to choose if plot of equation Jdot*qdot is equal to 0 during all simulation
               % (Jdot_qdotFlag == YES)> Outputs a plot of this equation
               % (Jdot_qdotFlag == NO)> No output of this plot 
               
%% Setting up simulation parameters

load('a4_Tsim.mat')
tf = Tsim;      % Final simulation time
h = 0.001;    % Number of samples within final vectors (times and states, sampling time)
opts = odeset('MaxStep',1e-3); % Max integration time step

%% Setting up model parameters & matrices
load("m_struc")
load("sm_struc")

%% Initial conditions
load('Otbot1_DataStates.mat')
% Initial conditions for the system
xs0 = x_out_data(1,:)';
clearvars x_out_data

%% Building the diferential equation
% Loading the data for the u input
load('Otbot1_Data.mat')     % Here we have the u's
load('a5_SampInstants.mat') % Here we have the sampling instants
load('a3_Tsaction.mat')     % Here we have the time of holding the action

dxdt= @(t,xs) xdot_otbot1_CTC_DF_IMU(t, xs, m, sm, FrictFlag, DistFlag, u_out_data, SampInst, Tsaction);

%% Simulationg with ode45
t = 0:h:tf;
[times,states]=ode45(dxdt,t,xs0,opts);

%% Compute KinEnergy in Every instant
if strcmp(KinEnFlag,'YES')
    ls = length(states);
    kinenvalues = zeros(ls,1);
   for i=1:ls
       kinenvalues(i,1) = sm.Texpr(states(i,3),states(i,9),states(i,6),states(i,11),states(i,12),states(i,10),states(i,7),states(i,8));
   end
end

%% Compute action torques in every instant

if strcmp(UvecFlag,'YES')
    u_vector = zeros(tf/h+1,4);
    ls = length(states);
    for i=1:ls
        u_f= u_function(times(i), u_out_data, SampInst, Tsaction);
        u_zoh = u_f;
        u_vector(i,1:4) = [u_zoh',times(i)];
    end
else
    disp('Plots of the evolution of actions u will not be displayed')
end

%% Compute accelerations and IMU readings
ls = length(states);
qdotdot_out = zeros(6,ls);
IMUxy = zeros(2,ls);
for i=1:ls
    Ms_out = sm.Mmatrix(states(i,3),states(i,6));
    Cs_out = sm.Cmatrix(states(i,3), states(i,9), states(i,6), states(i,12));
    Js_out = sm.Jmatrix(states(i,3), states(i,6));
    Es_out = sm.Ematrix;
    Jdots_out = sm.Jdotmatrix(states(i,3), states(i,9), states(i,6), states(i,12));

    % Setting Distrurbances
    switch DistFlag
        case 'YES'
            D_vec_out = dist_funct(times(i));
        case 'NO'
            D_vec_out = zeros(6,1);
        otherwise
            D_vec_out = zeros(6,1);
    end

    switch FrictFlag
        case 'YES'
            tau_friction_out(1:3,1) = zeros(3,1);
            tau_friction_out(4,1) = -m.b_frict(1,1)*states(i,10); % tau_friction_right (only works if we are in the .m file here we must use states(i,10), and put it inside the for loop)
            tau_friction_out(5,1) = -m.b_frict(2,1)*states(i,11); % tau_friction_left
            tau_friction_out(6,1) = -m.b_frict(3,1)*states(i,12); % tau_friction_pivot
        case 'NO'
            tau_friction_out = zeros(6,1);
        otherwise
            disp('Warning: This FrictFlag does not exist. Sisem qdotdot_out will be computed with no fricction')
            tau_friction_out = zeros(6,1);
    end

    % Compute the accelerations to have IMU readings 
    qdotdot_out(:,i) = [eye(6),zeros(6,3)]*pinv([Ms_out,Js_out.'; Js_out, zeros(3,3)])*[Es_out*u_vector(i,1:3)' + tau_friction_out + D_vec_out - Cs_out*states(i,7:12)'; -Jdots_out*states(i,7:12)'];
    
    % Transform the accelerations into IMU Readings (old way, using time
    % derivatives)
    % IMUxy(:,i) = m.dIRzmat22(states(i,3),states(i,9))*[states(i,7);states(i,8)] + m.IRzmat22(states(i,3))*qdotdot_out(1:2,i);
    
    % Transform the accelerations into IMU Readings (new way, only a projection of the acceleration vector to another base)
    IMUxy(:,i) = m.IRzmat22(states(i,3))*qdotdot_out(1:2,i);
end
qdotdot_out = qdotdot_out';
IMUxy = IMUxy';

%% Plot results
load('Otbot1_DataStates.mat')
load('Otbot1_Data.mat')
load('a5_SampInstants.mat')
load('Otbot1_DataAccel.mat')

% Load the extra flags needed
Biasalphadot_Flag = Flags.Biasalphadot_Flag;
Biasxddot_Flag = Flags.Biasxddot_Flag;
Biasyddot_Flag = Flags.Biasyddot_Flag;

%---------------------------------------------------------------%
%-                    Plots Minimized                          -%
%---------------------------------------------------------------%

% figure;
% hold on
% plot(times,states(:,1))
% plot(SampInst,x_out_data(:,1),'.', 'MarkerSize',15)
% hold off
% xlabel('t(s)','Interpreter','latex'),ylabel('$x$(meters)','Interpreter','latex'),title('Configuration','Interpreter','latex')
% legend('Estimated Model','Real Model','Interpreter','latex')
% 
% figure;
% hold on
% plot(times,states(:,2))
% plot(SampInst,x_out_data(:,2),'.', 'MarkerSize',15)
% hold off
% xlabel('t(s)','Interpreter','latex'),ylabel('$y$(meters)','Interpreter','latex'),title('Configuration','Interpreter','latex')
% legend('Estimated Model','Real Model','Interpreter','latex')
% 
% figure;
% hold on
% plot(times,states(:,3))
% plot(SampInst,x_out_data(:,3),'.', 'MarkerSize',15)
% hold off
% xlabel('t(s)','Interpreter','latex'),ylabel('$\alpha$(radians)','Interpreter','latex'),title('Configuration','Interpreter','latex')
% legend('Estimated Model','Real Model','Interpreter','latex')
% 
% figure;
% hold on
% plot(times,states(:,4))
% plot(SampInst,x_out_data(:,4),'.', 'MarkerSize',15)
% hold off
% xlabel('t(s)','Interpreter','latex'),ylabel('$\varphi_r$(radians)','Interpreter','latex'),title('Configuration','Interpreter','latex')
% legend('Estimated Model','Real Model','Interpreter','latex')
% 
% figure;
% hold on
% plot(times,states(:,5))
% plot(SampInst,x_out_data(:,5),'.', 'MarkerSize',15)
% hold off
% xlabel('t(s)','Interpreter','latex'),ylabel('$\varphi_l$(radians)','Interpreter','latex'),title('Configuration','Interpreter','latex')
% legend('Estimated Model','Real Model','Interpreter','latex')
% 
% figure;
% hold on
% plot(times,states(:,6))
% plot(SampInst,x_out_data(:,6),'.', 'MarkerSize',15)
% hold off
% xlabel('t(s)','Interpreter','latex'),ylabel('$\varphi_p$(radians)','Interpreter','latex'),title('Configuration','Interpreter','latex')
% legend('Estimated Model','Real Model','Interpreter','latex')
% 
% figure;
% hold on
% plot(times,states(:,7))
% plot(SampInst,x_out_data(:,7),'.', 'MarkerSize',15)
% hold off
% xlabel('t(s)','Interpreter','latex'),ylabel('$\dot{x}$(meters/second)','Interpreter','latex'),title('Velocity','Interpreter','latex')
% legend('Estimated Model','Real Model','Interpreter','latex')
% 
% figure;
% hold on
% plot(times,states(:,8))
% plot(SampInst,x_out_data(:,8),'.', 'MarkerSize',15)
% hold off
% xlabel('t(s)','Interpreter','latex'),ylabel('$\dot{y}$(meters/second)','Interpreter','latex'),title('Velocity','Interpreter','latex')
% legend('Estimated Model','Real Model','Interpreter','latex')
% 
% switch Biasalphadot_Flag
%     case 'NO'
%         figure;
%         hold on
%         plot(times,states(:,9))
%         plot(SampInst,y_out_data(:,1),'.', 'MarkerSize',15)
%         hold off
%         xlabel('t(s)','Interpreter','latex'),ylabel('$\dot{\alpha}$(radians/second)','Interpreter','latex'),title('Velocity','Interpreter','latex')
%         legend('Estimated Model','Sampled Data','Interpreter','latex')
%     case 'YES'
%         figure;
%         hold on
%         plot(times,states(:,9))
%         plot(SampInst,y_out_data(:,1),'.', 'MarkerSize',15)
%         plot(SampInst,x_out_data(:,9),'-.')
%         hold off
%         xlabel('t(s)','Interpreter','latex'),ylabel('$\dot{\alpha}$(radians/second)','Interpreter','latex'),title('Velocity','Interpreter','latex')
%         legend('Estimated Model','Sampled Data','Interpreter','latex')
%     otherwise
%         disp('WARNING: This Biasalphadot_Flag does not exist, omitting Rality plot')
%         figure;
%         hold on
%         plot(times,states(:,9))
%         plot(SampInst,y_out_data(:,1),'.', 'MarkerSize',15)
%         hold off
%         xlabel('t(s)','Interpreter','latex'),ylabel('$\dot{\alpha}$(radians/second)','Interpreter','latex'),title('Velocity','Interpreter','latex')
%         legend('Estimated Model','Sampled Data','Interpreter','latex')
% end
% 
% 
% figure;
% hold on
% plot(times,states(:,10))
% plot(SampInst,x_out_data(:,10),'.', 'MarkerSize',15)
% hold off
% xlabel('t(s)','Interpreter','latex'),ylabel('$\dot{\varphi}_r$(radians/second)','Interpreter','latex'),title('Velocity','Interpreter','latex')
% legend('Estimated Model','Real Model','Interpreter','latex')
% 
% figure;
% hold on
% plot(times,states(:,11))
% plot(SampInst,x_out_data(:,11),'.', 'MarkerSize',15)
% hold off
% xlabel('t(s)','Interpreter','latex'),ylabel('$\dot{\varphi}_l$(radians/second)','Interpreter','latex'),title('Velocity','Interpreter','latex')
% legend('Estimated Model','Real Model','Interpreter','latex')
% 
% figure;
% hold on
% plot(times,states(:,12))
% plot(SampInst,x_out_data(:,12),'.', 'MarkerSize',15)
% hold off
% xlabel('t(s)','Interpreter','latex'),ylabel('$\dot{\varphi}_p$(radians/second)','Interpreter','latex'),title('Velocity','Interpreter','latex')
% legend('Estimated Model','Real Model','Interpreter','latex')
% 
% % Plot accelerations x_ddot y_ddot  
% 
% switch Biasxddot_Flag
%     case 'NO'
%         figure;
%         hold on
%         plot(times,IMUxy(:,1))
%         plot(SampInst,y_out_data(:,2),'.', 'MarkerSize',15)
%         hold off
%         xlabel('t(s)','Interpreter','latex'),ylabel('$\ddot{x}$(meters/$s^2$)','Interpreter','latex'),title('Acceleration','Interpreter','latex')
%         legend('Estimated Model','Sampled Data','Interpreter','latex')
%     case 'YES'
%         figure;
%         hold on
%         plot(times,IMUxy(:,1))
%         plot(SampInst,y_out_data(:,2),'.', 'MarkerSize',15)
%         plot(SampInst,qddot_real_data(:,1),'-.')
%         hold off
%         xlabel('t(s)','Interpreter','latex'),ylabel('$\ddot{x}$(meters/$s^2$)','Interpreter','latex'),title('Acceleration','Interpreter','latex')
%         legend('Estimated Model','Sampled Data','Interpreter','latex')
%     otherwise
%         disp('WARNING: This Biasxddot_Flag does not exist, omitting Rality plot')
%         figure;
%         hold on
%         plot(times,IMUxy(:,1))
%         plot(SampInst,y_out_data(:,2),'.', 'MarkerSize',15)
%         hold off
%         xlabel('t(s)','Interpreter','latex'),ylabel('$\ddot{x}$(meters/$s^2$)','Interpreter','latex'),title('Acceleration','Interpreter','latex')
%         legend('Estimated Model','Sampled Data','Interpreter','latex')
% end
% 
% switch Biasyddot_Flag
%     case 'NO'
%         figure;
%         hold on
%         plot(times,IMUxy(:,2))
%         plot(SampInst,y_out_data(:,3),'.', 'MarkerSize',15)
%         hold off
%         xlabel('t(s)','Interpreter','latex'),ylabel('$\ddot{y}$(meters/$s^2$)','Interpreter','latex'),title('Acceleration','Interpreter','latex')
%         legend('Estimated Model','Sampled Data','Interpreter','latex')
%     case 'YES'
%         figure;
%         hold on
%         plot(times,IMUxy(:,2))
%         plot(SampInst,y_out_data(:,3),'.', 'MarkerSize',15)
%         plot(SampInst,qddot_real_data(:,2),'-.')
%         hold off
%         xlabel('t(s)','Interpreter','latex'),ylabel('$\ddot{y}$(meters/$s^2$)','Interpreter','latex'),title('Acceleration','Interpreter','latex')
%         legend('Estimated Model','Sampled Data','Interpreter','latex')
%     otherwise
%         disp('WARNING: This Biasyddot_Flag does not exist, omitting Rality plot')
%         figure;
%         hold on
%         plot(times,IMUxy(:,2))
%         plot(SampInst,y_out_data(:,3),'.', 'MarkerSize',15)
%         hold off
%         xlabel('t(s)','Interpreter','latex'),ylabel('$\ddot{y}$(meters/$s^2$)','Interpreter','latex'),title('Acceleration','Interpreter','latex')
%         legend('Estimated Model','Sampled Data','Interpreter','latex')
% end
% 
% % Plot the action vectors
% if strcmp(UvecFlag,'YES')
%     figure;
%     hold on
%     plot(times,u_vector(:,1:2))
%     plot(SampInst,u_out_data(:,1:2),'.', 'MarkerSize',15)
%     hold off
%     xlabel('t(s)','Interpreter','latex'),ylabel('$\tau$ of wheels(N$\cdot$m)','Interpreter','latex'),title('Action Torques','Interpreter','latex')
%     legend('$\tau$ right EM','$\tau$ left EM', '$\tau$ right SD', '$\tau$ left SD','Interpreter','latex')
% 
%     figure;
%     hold on
%     plot(times,u_vector(:,3))
%     plot(SampInst,u_out_data(:,3),'.', 'MarkerSize',15)
%     hold off
%     xlabel('t(s)','Interpreter','latex'),ylabel('$\tau$ pivot(N$\cdot$m)','Interpreter','latex'),title('Action Torques','Interpreter','latex')
%     legend('Estimated Model','Sampled Data','Interpreter','latex')
% end
% 
% if strcmp(KinEnFlag,'YES')
%     figure;
%     plot(times,kinenvalues)
%     xlabel('t(s)','Interpreter','latex'),ylabel('Kinetic Energy (Joules)','Interpreter','latex'), title('System Energy','Interpreter','latex')
% end


%---------------------------------------------------------------%
%-                    Plots Maximized                          -%
%---------------------------------------------------------------%
figure('WindowState','maximized');
hold on
plot(times,states(:,1))
plot(SampInst,x_out_data(:,1),'.', 'MarkerSize',15)
hold off
xlabel('t(s)','Interpreter','latex','FontSize', 17),ylabel('$x$(meters)','Interpreter','latex','FontSize', 17),title('Configuration','Interpreter','latex','FontSize', 17)
legend('Estimated Model','Real Model','Interpreter','latex','FontSize', 17)
grid on
box on


figure('WindowState','maximized');
hold on
plot(times,states(:,2))
plot(SampInst,x_out_data(:,2),'.', 'MarkerSize',15)
hold off
xlabel('t(s)','Interpreter','latex','FontSize', 17),ylabel('$y$(meters)','Interpreter','latex','FontSize', 17),title('Configuration','Interpreter','latex','FontSize', 17)
legend('Estimated Model','Real Model','Interpreter','latex','FontSize', 17)
grid on
box on

figure('WindowState','maximized');
hold on
plot(times,states(:,3))
plot(SampInst,x_out_data(:,3),'.', 'MarkerSize',15)
hold off
xlabel('t(s)','Interpreter','latex','FontSize', 17),ylabel('$\alpha$(radians)','Interpreter','latex','FontSize', 17),title('Configuration','Interpreter','latex','FontSize', 17)
legend('Estimated Model','Real Model','Interpreter','latex','FontSize', 17)
grid on
box on

figure('WindowState','maximized');
hold on
plot(times,states(:,4))
plot(SampInst,x_out_data(:,4),'.', 'MarkerSize',15)
hold off
xlabel('t(s)','Interpreter','latex','FontSize', 17),ylabel('$\varphi_r$(radians)','Interpreter','latex','FontSize', 17),title('Configuration','Interpreter','latex','FontSize', 17)
legend('Estimated Model','Real Model','Interpreter','latex','FontSize', 17)
grid on
box on

figure('WindowState','maximized');
hold on
plot(times,states(:,5))
plot(SampInst,x_out_data(:,5),'.', 'MarkerSize',15)
hold off
xlabel('t(s)','Interpreter','latex','FontSize', 17),ylabel('$\varphi_l$(radians)','Interpreter','latex','FontSize', 17),title('Configuration','Interpreter','latex','FontSize', 17)
legend('Estimated Model','Real Model','Interpreter','latex','FontSize', 17)
grid on
box on

figure('WindowState','maximized');
hold on
plot(times,states(:,6))
plot(SampInst,x_out_data(:,6),'.', 'MarkerSize',15)
hold off
xlabel('t(s)','Interpreter','latex','FontSize', 17),ylabel('$\varphi_p$(radians)','Interpreter','latex','FontSize', 17),title('Configuration','Interpreter','latex','FontSize', 17)
legend('Estimated Model','Real Model','Interpreter','latex','FontSize', 17)
grid on
box on

figure('WindowState','maximized');
hold on
plot(times,states(:,7))
plot(SampInst,x_out_data(:,7),'.', 'MarkerSize',15)
hold off
xlabel('t(s)','Interpreter','latex','FontSize', 17),ylabel('$\dot{x}$(meters/second)','Interpreter','latex','FontSize', 17),title('Velocity','Interpreter','latex','FontSize', 17)
legend('Estimated Model','Real Model','Interpreter','latex','FontSize', 17)
grid on
box on

figure('WindowState','maximized');
hold on
plot(times,states(:,8))
plot(SampInst,x_out_data(:,8),'.', 'MarkerSize',15)
hold off
xlabel('t(s)','Interpreter','latex','FontSize', 17),ylabel('$\dot{y}$(meters/second)','Interpreter','latex','FontSize', 17),title('Velocity','Interpreter','latex','FontSize', 17)
legend('Estimated Model','Real Model','Interpreter','latex','FontSize', 17)
grid on
box on

switch Biasalphadot_Flag
    case 'NO'
        figure('WindowState','maximized');
        hold on
        plot(times,states(:,9))
        plot(SampInst,y_out_data(:,1),'.', 'MarkerSize',15)
        hold off
        xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\dot{\alpha}$ [rad/s]','Interpreter','latex','FontSize', 17),title('Fitting plot for the angular velocity $\dot{\alpha}$','Interpreter','latex','FontSize', 17)
        legend('Estimated Model','Sampled Data','Interpreter','latex','FontSize', 17)
        grid on
        box on
    case 'YES'
        figure('WindowState','maximized');
        hold on
        plot(times,states(:,9))
        plot(SampInst,y_out_data(:,1),'.', 'MarkerSize',15)
        plot(SampInst,x_out_data(:,9),'-.')
        hold off
        xlabel('t(s)','Interpreter','latex','FontSize', 17),ylabel('$\dot{\alpha}$(radians/second)','Interpreter','latex','FontSize', 17),title('Velocity','Interpreter','latex','FontSize', 17)
        legend('Estimated Model','Sampled Data','Interpreter','latex','FontSize', 17)
        grid on
        box on
    otherwise
        disp('WARNING: This Biasalphadot_Flag does not exist, omitting Rality plot')
        figure('WindowState','maximized');
        hold on
        plot(times,states(:,9))
        plot(SampInst,y_out_data(:,1),'.', 'MarkerSize',15)
        hold off
        xlabel('t(s)','Interpreter','latex','FontSize', 17),ylabel('$\dot{\alpha}$(radians/second)','Interpreter','latex','FontSize', 17),title('Velocity','Interpreter','latex','FontSize', 17)
        legend('Estimated Model','Sampled Data','Interpreter','latex','FontSize', 17)
        grid on
        box on
end


figure('WindowState','maximized');
hold on
plot(times,states(:,10))
plot(SampInst,x_out_data(:,10),'.', 'MarkerSize',15)
hold off
xlabel('t(s)','Interpreter','latex','FontSize', 17),ylabel('$\dot{\varphi}_r$(radians/second)','Interpreter','latex','FontSize', 17),title('Velocity','Interpreter','latex','FontSize', 17)
legend('Estimated Model','Real Model','Interpreter','latex','FontSize', 17)
grid on
box on

figure('WindowState','maximized');
hold on
plot(times,states(:,11))
plot(SampInst,x_out_data(:,11),'.', 'MarkerSize',15)
hold off
xlabel('t(s)','Interpreter','latex','FontSize', 17),ylabel('$\dot{\varphi}_l$(radians/second)','Interpreter','latex','FontSize', 17),title('Velocity','Interpreter','latex','FontSize', 17)
legend('Estimated Model','Real Model','Interpreter','latex','FontSize', 17)
grid on
box on

figure('WindowState','maximized');
hold on
plot(times,states(:,12))
plot(SampInst,x_out_data(:,12),'.', 'MarkerSize',15)
hold off
xlabel('t(s)','Interpreter','latex','FontSize', 17),ylabel('$\dot{\varphi}_p$(radians/second)','Interpreter','latex','FontSize', 17),title('Velocity','Interpreter','latex','FontSize', 17)
legend('Estimated Model','Real Model','Interpreter','latex','FontSize', 17)
grid on
box on

% Plot accelerations x_ddot y_ddot  

switch Biasxddot_Flag
    case 'NO'
        figure('WindowState','maximized');
        hold on
        plot(times,IMUxy(:,1))
        plot(SampInst,y_out_data(:,2),'.', 'MarkerSize',15)
        hold off
        xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\ddot{x}$ [m/$\mathrm{s}^2$]','Interpreter','latex','FontSize', 17),title('Fitting plot for the acceleration $\ddot{x}$','Interpreter','latex','FontSize', 17)
        legend('Estimated Model','Sampled Data','Interpreter','latex','FontSize', 17)
        grid on
        box on
    case 'YES'
        figure('WindowState','maximized');
        hold on
        plot(times,IMUxy(:,1))
        plot(SampInst,y_out_data(:,2),'.', 'MarkerSize',15)
        plot(SampInst,qddot_real_data(:,1),'-.')
        hold off
        xlabel('t(s)','Interpreter','latex','FontSize', 17),ylabel('$\ddot{x}$(meters/$s^2$)','Interpreter','latex','FontSize', 17),title('Acceleration','Interpreter','latex','FontSize', 17)
        legend('Estimated Model','Sampled Data','Interpreter','latex','FontSize', 17)
        grid on
        box on
    otherwise
        disp('WARNING: This Biasxddot_Flag does not exist, omitting Rality plot')
        figure('WindowState','maximized');
        hold on
        plot(times,IMUxy(:,1))
        plot(SampInst,y_out_data(:,2),'.', 'MarkerSize',15)
        hold off
        xlabel('t(s)','Interpreter','latex','FontSize', 17),ylabel('$\ddot{x}$(meters/$s^2$)','Interpreter','latex','FontSize', 17),title('Acceleration','Interpreter','latex','FontSize', 17)
        legend('Estimated Model','Sampled Data','Interpreter','latex','FontSize', 17)
        grid on
        box on
end

switch Biasyddot_Flag
    case 'NO'
        figure('WindowState','maximized');
        hold on
        plot(times,IMUxy(:,2))
        plot(SampInst,y_out_data(:,3),'.', 'MarkerSize',15)
        hold off
        xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\ddot{y}$ [m/$\mathrm{s}^2$]','Interpreter','latex','FontSize', 17),title('Fitting plot for the acceleration $\ddot{y}$','Interpreter','latex','FontSize', 17)
        legend('Estimated Model','Sampled Data','Interpreter','latex','FontSize', 17)
        grid on
        box on
    case 'YES'
        figure('WindowState','maximized');
        hold on
        plot(times,IMUxy(:,2))
        plot(SampInst,y_out_data(:,3),'.', 'MarkerSize',15)
        plot(SampInst,qddot_real_data(:,2),'-.')
        hold off
        xlabel('t(s)','Interpreter','latex','FontSize', 17),ylabel('$\ddot{y}$(meters/$s^2$)','Interpreter','latex','FontSize', 17),title('Acceleration','Interpreter','latex','FontSize', 17)
        legend('Estimated Model','Sampled Data','Interpreter','latex','FontSize', 17)
        grid on
        box on
    otherwise
        disp('WARNING: This Biasyddot_Flag does not exist, omitting Rality plot')
        figure('WindowState','maximized');
        hold on
        plot(times,IMUxy(:,2))
        plot(SampInst,y_out_data(:,3),'.', 'MarkerSize',15)
        hold off
        xlabel('t(s)','Interpreter','latex','FontSize', 17),ylabel('$\ddot{y}$(meters/$s^2$)','Interpreter','latex','FontSize', 17),title('Acceleration','Interpreter','latex','FontSize', 17)
        legend('Estimated Model','Sampled Data','Interpreter','latex','FontSize', 17)
        grid on
        box on
end

% Plot the action vectors
if strcmp(UvecFlag,'YES')
    figure('WindowState','maximized');
    hold on
    plot(times,u_vector(:,1:2))
    plot(SampInst,u_out_data(:,1:2),'.', 'MarkerSize',15)
    hold off
    xlabel('t(s)','Interpreter','latex','FontSize', 17),ylabel('$\tau$ of wheels(N$\cdot$m)','Interpreter','latex','FontSize', 17),title('Action Torques','Interpreter','latex','FontSize', 17)
    legend('$\tau$ right EM','$\tau$ left EM', '$\tau$ right SD', '$\tau$ left SD','Interpreter','latex','FontSize', 17)
    grid on
    box on

    figure('WindowState','maximized');
    hold on
    plot(times,u_vector(:,3))
    plot(SampInst,u_out_data(:,3),'.', 'MarkerSize',15)
    hold off
    xlabel('t(s)','Interpreter','latex','FontSize', 17),ylabel('$\tau$ pivot(N$\cdot$m)','Interpreter','latex','FontSize', 17),title('Action Torques','Interpreter','latex','FontSize', 17)
    legend('Estimated Model','Sampled Data','Interpreter','latex','FontSize', 17)
    grid on
    box on
end

if strcmp(KinEnFlag,'YES')
    figure;
    plot(times,kinenvalues)
    xlabel('t(s)','Interpreter','latex'),ylabel('Kinetic Energy (Joules)','Interpreter','latex'), title('System Energy','Interpreter','latex')
end

%% Computing Jdot_qdot and the holonomic equation

if strcmp(Jdot_qdotFlag,'YES')
    ls = length(states);
    J_qvalues = zeros(ls,3);
    Hol_eq = zeros(ls,1);
    K_o = -(states(1,3) - m.r/(2*m.l_2)*states(1,4) + m.r/(2*m.l_2)*states(1,5) - states(1,6));
    for i=1:ls
        Jplot = sm.Jmatrix(states(i,3),states(i,6));
        qdotplot = states(i,7:12)';
        J_qvalues(i,:) = (Jplot*qdotplot)';
        Hol_eq(i) = states(i,3) - m.r/(2*m.l_2)*states(i,4) + m.r/(2*m.l_2)*states(i,5) - states(i,6) + K_o;
    end
    figure;
    hold on
    plot(times,J_qvalues)
    plot(times,Hol_eq)
    hold off
    xlabel('t(s)','Interpreter','latex'),ylabel('Result vector','Interpreter','latex'), title('Kinematic Error','Interpreter','latex')
    legend('$\mathbf{J} \cdot \dot{\mathbf{q}}_1$','$\mathbf{J} \cdot \dot{\mathbf{q}}_2$','$\mathbf{J} \cdot \dot{\mathbf{q}}_3$','Holeq','Interpreter','latex')
end

