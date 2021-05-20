clearvars
close all
clc
%% Flags
FlagVideo = 'YES'; % Flag to define if you want video output or not
               % (FlagVideo == YES)> simulation with VIDEO output
               % (FlagVideo == NO )> simulation without VIDEO output
               
u_Flag = 'CTE'; % Flag to define if you want a constnat torque for all simulation or not
               % (u_Flag == CTE)> simulation with constant torques
               % (u_Flag == VAR)> simulation with a torque function of time
               
TrackFlag = 'NO'; % Flag to define if you want to track the path of each center of mass
               % (TrackFlag == NO)> plots without any paths of center of
               % mass plotted
               
               % (TrackFlag == CB)> plots the path of the center of mass of
               % the chassis body
               
               % (TrackFlag == PF)> plots the path of the center of mass of
               % the platform body
               
               % (TrackFlag == BOTH)> plots the path of both centers of
               % mass

%% Setting up simulation parameters

tf = 10;      % Final simulation time
h = 0.001;    % Number of samples within final vectors (times and states)
opts = odeset('MaxStep',1e-3); % Max integration time step

%% Setting up model parameters & matrices
load("m_struc")
load("sm_struc")

%% Now we crate torques comand vector

u = [0.1;-0.1;0]; % All in Nm

%% Building the diferential equation
switch u_Flag
    case 'CTE'
        dxdt= @(t,xs) xdot_kin_otbot(t, xs, sm, u, u_Flag);
    case 'VAR'
        dxdt= @(t,xs) xdot_kin_otbot(t, xs, sm, u, u_Flag);
    otherwise
        FlagVideo = 'NO';
        disp('this flag does not exist')
end
   

% Initial conditions
xs0=zeros(12,1);

%% Simulationg with ode45
t = 0:h:tf;
[times,states]=ode45(dxdt,t,xs0,opts);

%% Plot results
figure;
plot(times,states(:,1))
xlabel('t(s)'),ylabel('x(meters)'),title('Configuration')

figure;
plot(times,states(:,2))
xlabel('t(s)'),ylabel('y(meters)'),title('Configuration')

figure;
plot(times,states(:,3))
xlabel('t(s)'),ylabel('alpha(radians)'),title('Configuration')

figure;
plot(times,states(:,4))
xlabel('t(s)'),ylabel('varphi_r(radians)'),title('Configuration')

figure;
plot(times,states(:,5))
xlabel('t(s)'),ylabel('varphi_l(radians)'),title('Configuration')

figure;
plot(times,states(:,6))
xlabel('t(s)'),ylabel('varphi_p(radians)'),title('Configuration')

figure;
plot(times,states(:,7))
xlabel('t(s)'),ylabel('xdot(meters/second)'),title('Velocity')

figure;
plot(times,states(:,8))
xlabel('t(s)'),ylabel('ydot(meters/second)'),title('Velocity')

figure;
plot(times,states(:,9))
xlabel('t(s)'),ylabel('alphadot(radians/second)'),title('Velocity')

figure;
plot(times,states(:,10))
xlabel('t(s)'),ylabel('varphidot_r(radians/second)'),title('Velocity')

figure;
plot(times,states(:,11))
xlabel('t(s)'),ylabel('varphidot_l(radians/second)'),title('Velocity')

figure;
plot(times,states(:,12))
xlabel('t(s)'),ylabel('varphidot_p(radians/second)'),title('Velocity')

%% Movie and video of the simulation
switch FlagVideo
    case 'YES'
        time=cputime;
        % Animation
        % h is the sampling time 
        % n is the scaling factor in order not to plot with the same step
        % than during the integration with ode45

        fs=30;
        n=round(1/(fs*h));

        % Set up the movie.
        writerObj = VideoWriter('otbot_sim','MPEG-4'); % Name it.
        %writerObj.FileFormat = 'mp4';
        writerObj.FrameRate = fs; % How many frames per second.
        open(writerObj);

        for i = 1:n:length(times)
                q.x = states(i,1);
                q.y = states(i,2);
                q.alpha = states(i,3);
                q.varphi_r = states(i,4);
                q.varphi_l = states(i,5);
                q.varphi_p = states(i,6);

                elapsed = cputime-time;
                if elapsed > 200
                    disp(elapsed);
                    disp('took too long to generate the video')
                    break
                else
                    draw_otbot(m,q)
                    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
                    writeVideo(writerObj, frame); 
                end

        end
        close(writerObj); % Saves the movie.
    case 'NO'
        disp('Simulation without video launched');
    otherwise
        disp('This flag input does not exist');
end
