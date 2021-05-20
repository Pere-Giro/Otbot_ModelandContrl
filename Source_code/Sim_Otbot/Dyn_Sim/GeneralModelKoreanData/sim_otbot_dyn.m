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
               
TrackFlag = 'PF'; % Flag to define if you want to track the path of each center of mass
               % (TrackFlag == NO)> plots without any paths of center of
               % mass plotted
               
               % (TrackFlag == CB)> plots the path of the center of mass of
               % the chassis body
               
               % (TrackFlag == PF)> plots the path of the center of mass of
               % the platform body
               
               % (TrackFlag == BOTH)> plots the path of both centers of
               % mass
               
KinEnFlag = 'YES'; % Flag to choose if Kinetic energic plot must be displayed or not
               % (KinEnFlag == YES)> Outputs a plot of kinetic energy
               % (KinEnFlag == NO )> No output plot of kinetic energy
%% Setting up simulation parameters

tf = 40;      % Final simulation time
h = 0.001;    % Number of samples within final vectors (times and states)
opts = odeset('MaxStep',1e-3); % Max integration time step

%% Setting up model parameters & matrices
load("m_struc")
load("sm_struc")

%% Initial conditions
p0 = [0,0,0]';
varphi0 = [0,0,0]';

varphidot0 = [5,1,0]';

% Need MFIK
MFIKmat = sm.MFIKmatrix(p0(3),varphi0(3));

% Computing p0
pdot0 = MFIKmat*varphidot0;

% Initial conditions for the system
xs0=[p0; varphi0; pdot0; varphidot0];

%% Now we crate torques comand vector

u = [0;0;0]; % All in Nm

%% Building the diferential equation
dxdt= @(t,xs) xdot_dyn_otbot(t, xs, sm, u, u_Flag);

%% Simulationg with ode45
t = 0:h:tf;
[times,states]=ode45(dxdt,t,xs0,opts);

%% Vectorize the path of centers of mass
if strcmp(TrackFlag,'NO')==0 
    
    % Setting this to adjust it to the video, here and down below they must
    % heve the same values, n & fs
    fs=30;
    n=round(1/(fs*h));
    
    % Set, ploting factor:
    Gpf = 20; % This plot factor is to scale if we want compute the com of each body every frame (Gpf = 1) or less frequently i.e every 2 frames (Gpf = 2...)
    n2 = Gpf*n;
    
    aux1=length(times);  
    GBpmat = zeros(3,round(aux1/n2));
    GPpmat = zeros(3,round(aux1/n2));
    jaux=1;

    for i = 1:n2:aux1
        q.x = states(i,1);
        q.y = states(i,2);
        q.alpha = states(i,3);
        q.varphi_r = states(i,4);
        q.varphi_l = states(i,5);
        q.varphi_p = states(i,6);

        [GBpi,GPpi] = c_o_m_bodies(m,q);
        GBpmat(:,jaux) = GBpi;
        GPpmat(:,jaux) = GPpi;
        jaux=jaux+1;
    end
end

%% Compute KinEnergy in Every instant
if strcmp(KinEnFlag,'YES')
    ls = length(states);
    kinenvalues = zeros(ls,1);
    
    kinenrot_b = zeros(ls,1);
    kinenrot_l = zeros(ls,1);
    kinenrot_r = zeros(ls,1);
    kinenrot_p = zeros(ls,1);
    kinenrot_sys = zeros(ls,1);
    kinenrot_chas = zeros(ls,1);
    
    kinentra_b = zeros(ls,1);
    kinentra_l = zeros(ls,1);
    kinentra_r = zeros(ls,1);
    kinentra_p = zeros(ls,1);
    kinentra_sys = zeros(ls,1);
    kinentra_chas = zeros(ls,1);
    
   for i=1:ls
       % Total kynetic energy
       kinenvalues(i,1) = sm.Texpr(states(i,3),states(i,9),states(i,6),states(i,11),states(i,12),states(i,10),states(i,7),states(i,8));
       
       % Rotation of the chassis body
       kinenrot_b(i,1) = sm.Trotbase(states(i,9),states(i,12));
       
       % Rotation of the right wheel
       kinenrot_r(i,1) = sm.Trot_rightwheel_expr(states(i,9),states(i,12),states(i,10));
       
       % Rotation of the left wheel
       kinenrot_l(i,1) = sm.Trot_leftwheel_expr(states(i,9),states(i,11),states(i,12));
       
       % Rotation of the platform
       kinenrot_p(i,1) = sm.Trot_platform_expr(states(i,9));
       
       % Total Rotational Energy of the Otbot
       kinenrot_sys(i,1) = kinenrot_b(i,1) + kinenrot_r(i,1) + kinenrot_l(i,1) + kinenrot_p(i,1);
       kinenrot_chas(i,1) = kinenrot_b(i,1) + kinenrot_r(i,1) + kinenrot_l(i,1);
       
       %%%%
       
       % Translation of the chassis body
       kinentra_b(i,1) = sm.Ttrabase(states(i,7),states(i,8));
       
       % Translation of the right wheel
       kinentra_r(i,1) = sm.Ttra_rightwheel_expr(states(i,3),states(i,9),states(i,6),states(i,12),states(i,7),states(i,8));
       
       % Translation of the left wheel
       kinentra_l(i,1) = sm.Ttra_leftwheel_expr(states(i,3),states(i,9),states(i,6),states(i,12),states(i,7),states(i,8));
       
       % Translation of the platform
       kinentra_p(i,1) = sm.Ttra_platform_expr(states(i,7),states(i,8));
       
       % Total Translational Energy of the Otbot
       kinentra_sys(i,1) = kinentra_b(i,1) + kinentra_r(i,1) + kinentra_l(i,1) + kinentra_p(i,1);
       kinentra_chas(i,1) = kinentra_b(i,1) + kinentra_r(i,1) + kinentra_l(i,1);
   
   end
end
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

if strcmp(KinEnFlag,'YES')
    figure;
    plot(times,kinenvalues)
    axis([0, tf, 0, 1])
    xlabel('t(s)'),ylabel('Total Kinetic Energy (Joules)'), title('System Energy')
    
    figure;
    plot(times,kinenrot_sys)
%     axis([0, tf, 0, 1])
    xlabel('t(s)'),ylabel('Rotational Energy (Joules)'), title('Rotational Kinetic Energy of Otbot')
    
    figure;
    plot(times,kinentra_sys)
%     axis([0, tf, 0, 1])
    xlabel('t(s)'),ylabel('Translational Energy (Joules)'), title('Translational Kinetic Energy of Otbot')
    
    figure;
    plot(times,kinenrot_chas)
%     axis([0, tf, 0, 1])
    xlabel('t(s)'),ylabel('Rotational Energy (Joules)'), title('Rotational Kinetic Energy of the Chassis')
    
    figure;
    plot(times,kinentra_chas)
%     axis([0, tf, 0, 1])
    xlabel('t(s)'),ylabel('Translational Energy (Joules)'), title('Translational Kinetic Energy of the Chassis')
    
    figure;
    plot(times,kinenrot_b)
%     axis([0, tf, 0, 1])
    xlabel('t(s)'),ylabel('Rotational Energy (Joules)'), title('Chassis Body (No wheels) Rotational Energy')
    
    figure;
    plot(times,kinenrot_r)
%     axis([0, tf, 0, 1])
    xlabel('t(s)'),ylabel('Rotational Energy (Joules)'), title('Right Wheel Rotational Energy')
    
    figure;
    plot(times,kinenrot_l)
%     axis([0, tf, 0, 1])
    xlabel('t(s)'),ylabel('Rotational Energy (Joules)'), title('Left Wheel Rotational Energy')
    
    figure;
    plot(times,kinenrot_p)
%     axis([0, tf, 0, 1])
    xlabel('t(s)'),ylabel('Rotational Energy (Joules)'), title('Platform Rotational Energy')
    
    figure;
    plot(times,kinentra_b)
%     axis([0, tf, 0, 1])
    xlabel('t(s)'),ylabel('Translational Energy (Joules)'), title('Chassis Body (No Wheels) Translational Energy')
    
    figure;
    plot(times,kinentra_r)
%     axis([0, tf, 0, 1])
    xlabel('t(s)'),ylabel('Translational Energy (Joules)'), title('Right Wheel Translational Energy')
    
    figure;
    plot(times,kinentra_l)
%     axis([0, tf, 0, 1])
    xlabel('t(s)'),ylabel('Translational Energy (Joules)'), title('Left Wheel Translational Energy')
    
    figure;
    plot(times,kinentra_p)
%     axis([0, tf, 0, 1])
    xlabel('t(s)'),ylabel('Translational Energy (Joules)'), title('Platform Translational Energy')
end

%% Movie and video of the simulation
switch FlagVideo
    case 'YES'
        switch TrackFlag
            case 'NO'
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
                        if elapsed > 1000
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
                
            case 'BOTH'
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
                jaux = 1;
                jaux2 = 1;
                for i = 1:n:length(times)
                        q.x = states(i,1);
                        q.y = states(i,2);
                        q.alpha = states(i,3);
                        q.varphi_r = states(i,4);
                        q.varphi_l = states(i,5);
                        q.varphi_p = states(i,6);

                        elapsed = cputime-time;
                        if elapsed > 1000
                            disp(elapsed);
                            disp('took too long to generate the video')
                            break
                        else
                            draw_otbot(m,q)
                            if jaux2 == 1
                                GBpplot = GBpmat(:,1:jaux);
                                GPpplot = GPpmat(:,1:jaux);
                                multi_circles(GBpplot,m.l_2/6,'b')
                                hold on
                                multi_circles(GPpplot,m.l_2/6,'g')
                                jaux2 = jaux2 + 1;
                                jaux = jaux+1;
                            elseif jaux2>= Gpf
                                multi_circles(GBpplot,m.l_2/6,'b')
                                hold on
                                multi_circles(GPpplot,m.l_2/6,'g')
                                jaux2 = 1;
                            else
                                multi_circles(GBpplot,m.l_2/6,'b')
                                hold on
                                multi_circles(GPpplot,m.l_2/6,'g')
                                jaux2 = jaux2 + 1;
                            end
                            frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
                            writeVideo(writerObj, frame); 
                        end
                end
                close(writerObj); % Saves the movie.
            case 'CB'
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
                jaux = 1;
                jaux2 = 1;
                for i = 1:n:length(times)
                        q.x = states(i,1);
                        q.y = states(i,2);
                        q.alpha = states(i,3);
                        q.varphi_r = states(i,4);
                        q.varphi_l = states(i,5);
                        q.varphi_p = states(i,6);

                        elapsed = cputime-time;
                        if elapsed > 1000
                            disp(elapsed);
                            disp('took too long to generate the video')
                            break
                        else
                            draw_otbot(m,q)
                            if jaux2 == 1
                                GBpplot = GBpmat(:,1:jaux);
                                
                                multi_circles(GBpplot,m.l_2/6,'b')
                                
                                
                                jaux2 = jaux2 + 1;
                                jaux = jaux+1;
                            elseif jaux2>= Gpf
                                multi_circles(GBpplot,m.l_2/6,'b')
                                
                                
                                jaux2 = 1;
                            else
                                multi_circles(GBpplot,m.l_2/6,'b')
                                
                                jaux2 = jaux2 + 1;
                            end
                            frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
                            writeVideo(writerObj, frame); 
                        end
                end
                close(writerObj); % Saves the movie.
            case 'PF'
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
                jaux = 1;
                jaux2 = 1;
                for i = 1:n:length(times)
                        q.x = states(i,1);
                        q.y = states(i,2);
                        q.alpha = states(i,3);
                        q.varphi_r = states(i,4);
                        q.varphi_l = states(i,5);
                        q.varphi_p = states(i,6);

                        elapsed = cputime-time;
                        if elapsed > 1000
                            disp(elapsed);
                            disp('took too long to generate the video')
                            break
                        else
                            draw_otbot(m,q)
                            if jaux2 == 1
                                
                                GPpplot = GPpmat(:,1:jaux);
                               
                                multi_circles(GPpplot,m.l_2/6,'g')
                                jaux2 = jaux2 + 1;
                                jaux = jaux+1;
                            elseif jaux2>= Gpf
                                
                                multi_circles(GPpplot,m.l_2/6,'g')
                                jaux2 = 1;
                            else
                                
                                multi_circles(GPpplot,m.l_2/6,'g')
                                jaux2 = jaux2 + 1;
                            end
                            frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
                            writeVideo(writerObj, frame); 
                        end
                end
                close(writerObj); % Saves the movie.
            otherwise
                disp('This TrackFlag does not exist, omitting the video')
        end            
    case 'NO'
        disp('Simulation without video launched');
    otherwise
        disp('This FlagVideo input does not exist');
end
