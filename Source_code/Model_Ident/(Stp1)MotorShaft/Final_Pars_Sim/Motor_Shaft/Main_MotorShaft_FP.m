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

UvecFlag = Flags.UvecFlag; % Flag to choose if you want to compute and display the plots of action vector u
                % (UvecFlag == YES)> action vector u plots will be displayed
                % (UvecFlag == NO)> action vector u plots will not be displayed
               
%% Setting up simulation parameters

load('a4_Tsim.mat')
tf = Tsim;      % Final simulation time
h = 0.001;    % Number of samples within final vectors (times and states, sampling time)
opts = odeset('MaxStep',1e-3); % Max integration time step

%% Setting up model parameters & matrices
load("m_struc")

%% Initial conditions
load('MotorShaft_DataStates.mat')
% Initial conditions for the system
xs0 = x_out_data(1,:)';
clearvars x_out_data

%% Building the diferential equation
% Loading the data for the u input
load('MotorShaft_Data.mat')     % Here we have the u's
load('a5_SampInstants.mat')     % Here we have the sampling instants
load('a3_Tsaction.mat')         % Here we have the time of holding the action

dxdt= @(t,xs) xdot_motorshaft(t, xs, m, u_out_data, SampInst, Tsaction);

%% Simulationg with ode45
t = 0:h:tf;
[times,states]=ode45(dxdt,t,xs0,opts);

%% Compute action torques in every instant

if strcmp(UvecFlag,'YES')
    u_vector = zeros(tf/h+1,2);
    ls = length(states);
    for i=1:ls
        u_f= u_function(times(i), u_out_data, SampInst, Tsaction);
        u_zoh = u_f;
        u_vector(i,1:2) = [u_zoh',times(i)];
    end
else
    disp('Plots of the evolution of actions u will not be displayed')
end

%% Plot results
load('MotorShaft_DataStates.mat')
load('MotorShaft_Data.mat')
load('a5_SampInstants.mat')


%------------------------------------------------------------%
%- Plots tunned for the TFM (default colors dots over line) -%
%------------------------------------------------------------%

% figure('WindowState','maximized');
% hold on
% plot(times,states(:,1))
% plot(SampInst,x_out_data(:,1),'.','MarkerSize',15)
% hold off
% xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\varphi$ [rad]','Interpreter','latex','FontSize', 17),title('Configuration','Interpreter','latex','FontSize', 17)
% legend('Estimated Model','Real Model','Interpreter','latex','FontSize', 17)
% grid on
% box on
% 
% % Plot for the wheel [comment if pivot]
% % figure('WindowState','maximized');
% % hold on
% % plot(times,states(:,2))
% % plot(SampInst,y_out_data,'.','MarkerSize',15)
% % hold off
% % xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\dot{\varphi}$ [rad/s]','Interpreter','latex','FontSize', 17),title('Fitting plot for angular velocity of a wheel $\dot{\varphi}_i$','Interpreter','latex','FontSize', 17)
% % legend('Estimated Model','Sampled Data','Interpreter','latex','FontSize', 17)
% % grid on
% % box on
% 
% % Plot for the pivot [comment if wheel]
% figure('WindowState','maximized');
% hold on
% plot(times,states(:,2))
% plot(SampInst,y_out_data,'.','MarkerSize',15)
% hold off
% xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\dot{\varphi}$ [rad/s]','Interpreter','latex','FontSize', 17),title('Fitting plot for angular velocity of the pivot shaft $\dot{\varphi}_p$','Interpreter','latex','FontSize', 17)
% legend('Estimated Model','Sampled Data','Interpreter','latex','FontSize', 17)
% grid on
% box on
% 
% % Plot the action vectors
% if strcmp(UvecFlag,'YES')
%     figure('WindowState','maximized');
%     hold on
%     plot(times,u_vector(:,1))
%     plot(SampInst,u_out_data(:,1),'.','MarkerSize',15)
%     hold off
%     xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\tau$ [N$\cdot$m]','Interpreter','latex','FontSize', 17),title('Action Torque','Interpreter','latex','FontSize', 17)
%     legend('Estimated Model','Sampled Data','Interpreter','latex','FontSize', 17)
%     grid on
%     box on
% end




%------------------------------------------------------------%
%- Plots tunned for the TFM (red & green, line over dots)   -%
%------------------------------------------------------------%

figure('WindowState','maximized');
hold on
plot(SampInst,x_out_data(:,1),'r.','MarkerSize',15)
plot(times,states(:,1),'g')
hold off
xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\varphi$ [rad]','Interpreter','latex','FontSize', 17),title('Configuration','Interpreter','latex','FontSize', 17)
legend('Estimated Model','Real Model','Interpreter','latex','FontSize', 17)
grid on
box on

% Plot for the wheel [comment if pivot]
figure('WindowState','maximized');
hold on
plot(times,states(:,2),'g','LineWidth',2)
plot(SampInst,y_out_data,'r.','MarkerSize',15)
hold off
xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\dot{\varphi}_r$ [rad/s]','Interpreter','latex','FontSize', 17),title('Fitting plot for angular velocity of the right wheel $\dot{\varphi}_r$','Interpreter','latex','FontSize', 17)
legend('Estimated model','Sampled data','Interpreter','latex','FontSize', 17)
grid on
box on

% Plot for the pivot [comment if wheel]
% figure('WindowState','maximized');
% hold on
% plot(times,states(:,2),'g','LineWidth',2)
% plot(SampInst,y_out_data,'r.','MarkerSize',7)
% hold off
% xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\dot{\varphi}_p$ [rad/s]','Interpreter','latex','FontSize', 17),title('Fitting plot for angular velocity of the pivot shaft $\dot{\varphi}_p$','Interpreter','latex','FontSize', 17)
% legend('Estimated model','Sampled data','Interpreter','latex','FontSize', 17)
% grid on
% box on

% Plot the action vectors
if strcmp(UvecFlag,'YES')
    figure('WindowState','maximized');
    hold on
    plot(SampInst,u_out_data(:,1),'r.','MarkerSize',15)
    plot(times,u_vector(:,1),'g')
    hold off
    xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\tau$ [N$\cdot$m]','Interpreter','latex','FontSize', 17),title('Action Torque','Interpreter','latex','FontSize', 17)
    legend('Estimated Model','Sampled Data','Interpreter','latex','FontSize', 17)
    grid on
    box on
end


