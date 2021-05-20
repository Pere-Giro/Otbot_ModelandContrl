clearvars
close all
clc
%% Flags
FlagVideo = 'NO'; % Flag to define if you want video output or not
               % (Flag1 == YES)> simulation with VIDEO output
               % (Flag1 == NO )> simulation without VIDEO output

%% Setting up simulation parameters

tf = 15;      % Final simulation time
h = 0.01;    % Number of samples within final vectors (times and states)
% opts = odeset('MaxStep',1e-3); % Max integration time step

%% Setting up model parameters & matrices
load("m_struc")
load("sm_struc")

%% Now we crate the velocity comand vector

v = [1;1;1]; % All in rad/s

%% Building the diferential equation
dxdt= @(t,xs) qdot_kin_otbot(t, xs, sm, v);

% Initial conditions
xs0=zeros(6,1);

%% Simulationg with ode45
t = 0:h:tf;
[times,states]=ode45(dxdt,t,xs0);

%% Plot results
figure;
plot(times,states(:,1))
xlabel('t(s)'),ylabel('x(t)'),title('Configuration')

figure;
plot(times,states(:,2))
xlabel('t(s)'),ylabel('y(t)'),title('Configuration')

figure;
plot(times,states(:,3))
xlabel('t(s)'),ylabel('alpha(t)'),title('Configuration')

figure;
plot(times,states(:,4))
xlabel('t(s)'),ylabel('varphi_r(t)'),title('Configuration')

figure;
plot(times,states(:,5))
xlabel('t(s)'),ylabel('varphi_l(t)'),title('Configuration')

figure;
plot(times,states(:,6))
xlabel('t(s)'),ylabel('varphi_p(t)'),title('Configuration')

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
