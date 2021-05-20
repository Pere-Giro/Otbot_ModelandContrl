clearvars
close all
clc
%% WARNING
% This simulation uses allways a zoh to give the input action to the
% system, meaning that the actions are not updated at the same frequancy
% that the evolution of the system itself
%% Flags
               
Noise_Flag = 'YES'; % Flag to define if you want to set noise in the data sampling of alphadot
                                % (Noise_Flag == YES)> There will be noise added
                                % to sampling of varphidot
                                % (Noise_Flag == NO )> There will be no noise
                                % added to sampling of varphidot
                        
                       SENSOR.Stdev = 0.01;                  
                                     
u_Flag = 'CTE'; % Flag to define if you want a constnat torque for all simulation or not
               % (u_Flag == CTE)> simulation with constant torques
               % (u_Flag == VAR)> simulation with a torque function of time
                
                % Now we create initial torques comand vector
                u = 0; % All in Nm
                u_cte = 6; % This is the torque vector used if Flag is set to CTE
               
UvecFlag = 'YES'; % Flag to choose if you want to compute and display the plots of action vector u
                % (UvecFlag == YES)> action vector u plots will be displayed
                % (UvecFlag == NO)> action vector u plots will not be displayed
                           
%% Setting up simulation parameters

tf = 0.5;      % Final simulation time
h = 0.001;    % Number of samples within final vectors (times and states, sampling time)
opts = odeset('MaxStep',1e-3); % Max integration time step

%% Setting up model parameters & matrices
load("m_struc")

%% ZOH time variable
h2 = 0.01;    % This is the sampling time for having our reality data sampled
h3 = 0.01;    % This is the time interval during which torque actions will remain constant
global tstep
tstep.t = 0;  % Variable to have a time reference for constant actuation during steps initial value is set to the value of sampling time
tstep.h = h2; % Also give the value of sampling time for the data
tstep.haction = h3; % Time interval to hold the same action (zoh)
tstep.u = u;  % We also give the value of initial u vector this will keep updateing
%% Initial conditions

xs0 = zeros(2,1);

%% Building the diferential equation
dxdt= @(t,xs) xdot_motorshaft(t, xs, m, u, u_cte, u_Flag);

%% Simulationg with ode45
t = 0:h:tf;
[times,states]=ode45(dxdt,t,xs0,opts);

%% Compute action torques in every instant
if strcmp(UvecFlag,'YES')
    u_vector = zeros(tf/h+1,2);
    Int_Vecout = zeros(1,7); %[I_k-1(1,3), error_k-1(1,3), t_k-1(1)]
    % Resetting tstep
    tstep.t = 0;        % Variable to have a time reference for constant actuation during steps initial value is set to the value of sampling time
    tstep.h = h2;       % Also give the value of sampling time
    tstep.haction = h3; % Time interval to hold the same action (zoh)
    tstep.u = u;        % We also give the value of initial u vector this will keep updateing
    switch u_Flag
        case 'CTE'
            ls = length(states);
            for i=1:ls
                u_vector(i,1:2) = [u_cte',times(i)];
            end
        case 'VAR'
            ls = length(states);
            for i=1:ls
                u_f= u_function(times(i), u);
                u_zoh = zoh_function(times(i),u_f);
                u_vector(i,1:2) = [u_zoh',times(i)];
            end
        otherwise
            disp('Warning this u_Flag does not exits omiting torques camputation')
            disp('Computation with CTE torques will be launched')
            ls = length(states);
            for i=1:ls
                u_vector(i,1:2) = [u_cte',times(i)];
            end
    end
else
    disp('Plots of the evolution of actions u will not be displayed')
end

%% Plot results
figure;
plot(times,states(:,1))
xlabel('t(s)','Interpreter','latex'),ylabel('$\varphi$(rad)','Interpreter','latex'),title('Configuration','Interpreter','latex')

figure;
plot(times,states(:,2))
xlabel('t(s)','Interpreter','latex'),ylabel('$\dot{\varphi}$(rad/s)','Interpreter','latex'),title('Velocity','Interpreter','latex')

% Plot the action vectors
if strcmp(UvecFlag,'YES')
    figure;
    plot(times,u_vector(:,1))
    xlabel('t(s)','Interpreter','latex'),ylabel('$\tau$(N$\cdot$m)','Interpreter','latex'),title('Action Torque','Interpreter','latex')
end

%% Create output data
Upath = userpath;

ls = length(states);
if h2>h
    hts = uint64(h2/h);
    ls = length(states);
    
    % Build data output y including alpha_dot, x_ddot, y_ddot
    % First we add the proportional Bias aka Scaling factor error
    y_out_data1 = states(1:hts:ls,2);
    % Then we add the noise
    y_out_data = add_noise(Noise_Flag, y_out_data1, SENSOR);
    
    % Build action output 
    u_out_data = u_vector(1:hts:ls,1);
    
    % Save the real values of the parameters
    real_pars = realpfunc(m);

    % Save the values of the state & accelerations to compare it later
    x_out_data = states(1:hts:ls,:);
    
    savedir1 = strcat(Upath,'\Model_Identification\MotorShaft(Stp1)\');
    savedir2 = strcat(Upath,'\Model_Identification\MotorShaft(Stp1)\Final_Pars_Sim\Motor_Shaft\');
    
    save(strcat(savedir1,'MotorShaft_Data.mat'),"y_out_data","u_out_data")
    save(strcat(savedir2,'MotorShaft_Data.mat'),"y_out_data","u_out_data")
    save(strcat(savedir2,'MotorShaft_DataStates.mat'),"x_out_data")
    save(strcat(savedir1,'initial_states.mat'),"xs0")
    save(strcat(savedir1,'real_pars_data.mat'),"real_pars")
else
    disp('Warning h2 has to be always greater than h omiting output vectors')
end

%% Export the Flags
% Create a struc object to save all the flags and export them

Flags.Noise_Flag = Noise_Flag;             
Flags.u_Flag = u_Flag;                          
Flags.UvecFlag = UvecFlag;                         
Flags.u = u;
Flags.u_cte = u_cte;
Flags.Stdev = SENSOR.Stdev;

save(strcat(savedir2,'Flags_Data.mat'),"Flags")

