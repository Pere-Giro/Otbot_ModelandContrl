clearvars
close all
clc
%% WARNING
% This simulation uses allways a zoh to give the input action to the
% system, meaning that the actions are not updated at the same frequancy
% that the evolution of the system itself
%% Flags
FlagVideo = 'NO'; % Flag to define if you want video output or not
               % (FlagVideo == YES)> simulation with VIDEO output
               % (FlagVideo == NO )> simulation without VIDEO output
               
Noisealphadot_Flag = 'YES'; % Flag to define if you want to set noise in the data sampling of alphadot
                       % (Biasalphadot_Flag == YES)> There will be noise added
                       % to sampling of alphadot
                       % (Biasalphadot_Flag == NO )> There will be no noise
                       % added to sampling of alphadot
                       
Nopt_Flag1 = 'DATA'; % Flag to decide if the stdev of the Noise of alphadot is computed using
                     % SR & ND or is set to a value directly
                        % (Nopt_Flag1 == DATA )> StDev is computed using the data
                        % (Nopt_Falg1 == VALUE)> StDev is set to the value
                        % directly
                                            
                       % Set the data to compute stdev
                       IMU.ASR = 0.01^-1; % [Hz]
                       IMU.NDgy = 0.14e-3*9.81; % mg/sqrt(Hz)
                       % Set the value of Stdev directly
                       IMU.StdevA = 0.05; 
               
Noisexddot_Flag = 'YES'; % Flag to define if you want to set noise in the data sampling of x_ddot
                       % (Biasxddot_Flag == YES)> There will be noise added
                       % to sampling of x_ddot
                       % (Biasxddot_Flag == NO )> There will be no noise
                       % added to sampling of x_ddot
                       
Nopt_Flag2 = 'DATA'; % Flag to decide if the stdev of the Noise of x_ddot is computed using
                     % SR & ND or is set to a value directly
                        % (Nopt_Flag2 == DATA )> StDev is computed using the data
                        % (Nopt_Falg2 == VALUE)> StDev is set to the value
                        % directly
                                            
                       % Set the data to compute stdev
                       IMU.XSR = 0.01^-1; % [Hz]
                       IMU.NDaccX = 0.14e-3*9.81; % mg/sqrt(Hz)
                       % Set the value of Stdev directly
                       IMU.StdevX = 0.05;
                       
Noiseyddot_Flag = 'YES'; % Flag to define if you want to set noise in the data sampling of y_ddot
                       % (Biasyddot_Flag == YES)> There will be noise added
                       % to sampling of y_ddot
                       % (Biasyddot_Flag == NO )> There will be no noise
                       % added to sampling of y_ddot
                       
Nopt_Flag3 = 'DATA'; % Flag to decide if the stdev of the Noise of y_ddot is computed using
                     % SR & ND or is set to a value directly
                        % (Nopt_Flag3 == DATA )> StDev is computed using the data
                        % (Nopt_Falg3 == VALUE)> StDev is set to the value
                        % directly
                                            
                       % Set the data to compute stdev
                       IMU.YSR = 0.01^-1; % [Hz]
                       IMU.NDaccY = 0.14e-3*9.81; % mg/sqrt(Hz)
                       % Set the value of Stdev directly
                       IMU.StdevY = 0.05;  
                                     
u_Flag = 'CTE'; % Flag to define if you want a constnat torque for all simulation or not
               % (u_Flag == CTE)> simulation with constant torques
               % (u_Flag == VAR)> simulation with a torque function of time
               % (u_Flag == CTC_LQR)> simulation with a CTC controler using
               % LQR
               % (u_Flag == CTC_PP)> simulation with a CTC controler using
               % pole placement method
                
                % Now we create initial torques comand vector
                u = [0;0;0]; % All in Nm
                u_cte = [6;-10;6]; % This is the torque vector used if Flag is set to CTE
                
u_IntFlag = 'NO'; % Flag to define if you want a controler including an Integrative component or not
               % (u_IntFlag == YES) the CTC controller tuned via PP will
               % include an integrative term.
               % (u_IntFlag == NO) the CTC controller tuned via PP will not
               % include an integrative term.
               
goal_Flag = 'FIX'; % Flag to define the desired trajectory for the robot, only works if u is set to LQR or PP
               % (goal_Flag == FIX)> goal trajectory is a fixed point
               % (goal_Flag == CIRCLE)> goal trajectory is a circle
               % (goal_Flag == POLILINE)> goal trajectory is a poliline
               
               % Set a fixed point in case goal_Flag = FIX
               Fix_point = [0.2,0.2,deg2rad(10),0,0,0]'; % [x,y,alpha,varphi_r,varphi_l,varphi_p] velocities are 0
               
FrictFlag = 'NO'; % Flag to define if we want to add friction in our model or not
               % (FrictFlag == YES)> Our model will have freictions
               % (FrictFlag == NO)> Our model will not contain frictions
               
DistFlag = 'NO'; % Flag to define if we want to add non modelled disturbances during the simulation
               % (DistFlag == YES)> Non modelled disturbances will be added 
               % (DistFlag == NO)> Non modelled disturbances will NOT be added
               
TrackFlag = 'PF'; % Flag to define if you want to track the path of each center of mass
               % (TrackFlag == NO)> plots without any paths of center of
               % mass plotted              
               % (TrackFlag == CB)> plots the path of the center of mass of
               % the chassis body               
               % (TrackFlag == PF)> plots the path of the center of mass of
               % the platform body              
               % (TrackFlag == BOTH)> plots the path of both centers of
               % mass
               
TargetFlag = 'NO'; % Flag to choose if you want to see the desired goal in the video simulation
                    % (Only avaliable if u_Flag is CTC_LQR or CTC_PP)
                % (TargetFlag == YES)> A red target will be displayed in the
                % simulation
                % (TargetFlag == NO)> There won't be a target in the video
                % simulation
               
UvecFlag = 'YES'; % Flag to choose if you want to compute and display the plots of action vector u
                % (UvecFlag == YES)> action vector u plots will be displayed
                % (UvecFlag == NO)> action vector u plots will not be displayed
                           
KinEnFlag = 'NO'; % Flag to choose if Kinetic energic plot must be displayed or not
               % (KinEnFlag == YES)> Outputs a plot of kinetic energy
               % (KinEnFlag == NO )> No output plot of kinetic energy
               
Jdot_qdotFlag = 'YES'; % Flag to choose if plot of equation Jdot*qdot is equal to 0 during all simulation
               % (Jdot_qdotFlag == YES)> Outputs a plot of this equation
               % (Jdot_qdotFlag == NO)> No output of this plot 
               
ErrorFlag = 'NO'; % Flag to choose if error plots are displayed or not
               % (ErrorFlag == YES)> Error plots will be displayed
               % (ErrorFlag == NO)> Error plots will not be displayed
               
Biasalphadot_Flag = 'NO'; % Flag to define if you want to set a bias in the data sampling of alphadot
                       % (Biasalphadot_Flag == YES)> There will be a bias added
                       % to sampling of alphadot
                       % (Biasalphadot_Flag == NO )> There will be no bias
                       % added to sampling of alphadot
                       
                       % Set the desired bias 
                       IMU.Abias = 1.5;
               
Biasxddot_Flag = 'NO'; % Flag to define if you want to set a bias in the data sampling of x_ddot
                       % (Biasxddot_Flag == YES)> There will be a bias added
                       % to sampling of x_ddot
                       % (Biasxddot_Flag == NO )> There will be no bias
                       % added to sampling of x_ddot
                       
                       % Set the desired bias 
                       IMU.Xbias = 2;
                       
Biasyddot_Flag = 'NO'; % Flag to define if you want to set a bias in the data sampling of y_ddot
                       % (Biasyddot_Flag == YES)> There will be a bias added
                       % to sampling of y_ddot
                       % (Biasyddot_Flag == NO )> There will be no bias
                       % added to sampling of y_ddot
                       
                       % Set the desired bias 
                       IMU.Ybias = 3;
               
               
%% Setting up simulation parameters

tf = 3;      % Final simulation time
h = 0.001;    % Number of samples within final vectors (times and states, sampling time)
opts = odeset('MaxStep',1e-3); % Max integration time step

%% Setting up model parameters & matrices
load("m_struc")
load("sm_struc")

%% ZOH time variable
h2 = 0.01;    % This is the sampling time for having our reality data sampled
h3 = 0.01;    % This is the time interval during which torque actions will remain constant
global tstep
tstep.t = 0;  % Variable to have a time reference for constant actuation during steps initial value is set to the value of sampling time
tstep.h = h2; % Also give the value of sampling time for the data
tstep.haction = h3; % Time interval to hold the same action (zoh)
tstep.u = u;  % We also give the value of initial u vector this will keep updateing
%% Initial conditions
p0 = [0, 0, 0]';
varphi0 = [0, 0, 0]';

pdot0 = zeros(3,1);

% Initial condicions for the DirectCircle [OpenLoop Circle] (Makeing otbot follow a circular trajectory without controler only with u(t)
% p0 = [1, m.l_1, 0]';
% varphi0 = [0, 0, -pi/2]';
% 
% pdot0 = [-m.l_1, 1, 0]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Need MIIK
MIIKmat = sm.MIIKmatrix(p0(3),varphi0(3));
% Computing p0
varphidot0 = MIIKmat*pdot0;

% Initial conditions for the system
xs0=[p0; varphi0; pdot0; varphidot0];

%% Desired trajectory
if strcmp(u_Flag,'CTC_LQR') || strcmp(u_Flag,'CTC_PP')
    if strcmp(goal_Flag,'FIX')
        cp.pss = Fix_point;
    else
        cp.pss = zeros(6,1);
    end
end
% cp struct object stands for control parameters inside contains the desires
% goal fix point and the gain matrix

%% Creating gain matrices for the CTC 

% Our linear system matrices are
Actc = [zeros(3,3),eye(3,3); zeros(3,3), zeros(3,3)];
Bctc = [zeros(3,3);eye(3,3)];

if strcmp(u_Flag,'CTC_LQR')
    Qctc = [3e5*eye(3,3),zeros(3,3);
            zeros(3,3), 5*eye(3,3)];
    Rctc = 5*eye(3,3);
    cp.K = lqr(Actc,Bctc,Qctc,Rctc);
    
elseif strcmp(u_Flag,'CTC_PP')
    if strcmp(u_IntFlag,'NO')
%         EIG_set = [-1,-1.1,-1.2, -1.3, -1.4, -1.5];
        EIG_set = [-0.8,-8,-0.8, -8, -0.8, -8]; % This is the good set of poles in order to have a 5 second time of establishment
        cp.K = place(Actc,Bctc,EIG_set);
    elseif strcmp(u_IntFlag,'YES')
        Actc = [zeros(3,3),eye(3,3), zeros(3,3);
                zeros(3,3), zeros(3,3), eye(3,3);
                zeros(3,9)];
        Bctc = [zeros(6,3);eye(3,3)];
        EIG_set = [-0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1];
        cp.K = place(Actc,Bctc,EIG_set);   
    else
        disp('This u_IntFlag does not exist, tuning CTC as standart PD')
        EIG_set = [-1,-1.1,-1.2, -1.3, -1.4, -1.5];
        cp.K = place(Actc,Bctc,EIG_set);
    end
else
    cp.K =0; % This is to avoid getting Matlab error due to the fact thet cp does not exist
end

%% Setting-up the Int_Vec (Integrative Vector)
% This part is designed to create a vector that will contain the previous
% values of the Integral and it's time intant. We need this, beacuse we
% will be approximating the integral using the trapezoidal rule.

Int_Vec = zeros(1,7); %[I_k-1(1,3), error_k-1(1,3), t_k-1(1)]

%% Building the diferential equation
dxdt= @(t,xs) xdot_otbot1_CTC_DF_IMU(t, xs, m, sm, u, u_cte, u_Flag, cp, goal_Flag, FrictFlag, DistFlag, u_IntFlag, Int_Vec);

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
    Gpf = 16; % This plot factor is to scale if we want compute the com of each body every frame (Gpf = 1) or less frequently i.e every 2 frames (Gpf = 2...)
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
   for i=1:ls
       kinenvalues(i,1) = sm.Texpr(states(i,3),states(i,9),states(i,6),states(i,11),states(i,12),states(i,10),states(i,7),states(i,8));
   end
end

%% Compute the error for every instant
uswitch = strcmp(u_Flag,'CTC_LQR') || strcmp(u_Flag,'CTC_PP');

if strcmp(ErrorFlag,'YES') && uswitch
    lst = length(states);
    errorstates = zeros(lst,6);
    switch goal_Flag
        case 'FIX'
            for i = 1:lst
                errorstates(i,1:3) = cp.pss(1:3)' - states(i,1:3);
                errorstates(i,4:6) = cp.pss(4:6)' - states(i,7:9);
            end
        case 'CIRCLE'
            for i=1:lst
                [Xc,Yc,Alphac,Xdotc,Ydotc,Alphadotc,Xddotc,Yddotc,Alphaddotc] = TD_circle_of_time_XYA(times(i));
                errorstates(i,1:3) = [Xc,Yc,Alphac] - states(i,1:3);
                errorstates(i,4:6) = [Xdotc,Ydotc,Alphadotc] - states(i,7:9);
            end
        case 'POLILINE'
            for i=1:lst
                [Xc,Yc,Alphac,Xdotc,Ydotc,Alphadotc,Xddotc,Yddotc,Alphaddotc] = TD_polyline_of_time(times(i));
                errorstates(i,1:3) = [Xc,Yc,Alphac] - states(i,1:3);
                errorstates(i,4:6) = [Xdotc,Ydotc,Alphadotc] - states(i,7:9);
            end
            
    end
end

%% Compute desired GOAL (Target) in Every instant

uswitch = strcmp(u_Flag,'CTC_LQR') || strcmp(u_Flag,'CTC_PP');

if strcmp(TargetFlag,'YES') && strcmp(FlagVideo,'YES') && uswitch
    TEFlag = 0;
    lst = length(states);
    switch goal_Flag
        case 'FIX'
            targetpos = repmat(cp.pss',[lst,1]);
            
        case 'CIRCLE'
            targetpos = zeros(lst,3);
            for i=1:lst
                [Xc,Yc,Alphac,aux1,aux2,aux3,aux4,aux5,aux6] = TD_circle_of_time_XYA(times(i));
                targetpos(i,:) = [Xc,Yc,Alphac];
            end
        case 'POLILINE'
            targetpos = zeros(lst,3);
            for i=1:lst
                [Xc,Yc,Alphac,aux1,aux2,aux3,aux4,aux5,aux6] = TD_polyline_of_time(times(i));
                targetpos(i,:) = [Xc,Yc,Alphac];
            end
        otherwise
            disp('This goal_Flag does not exist cant compute the Target')
            TEFlag = 1;
    end                      
else 
    TEFlag = 1;
    disp('The combination of flags used is not appropriate to compute Target')
    disp('If you want, you can change the combination of flags in order to display the desired results')
    disp('Target location will not be computed')
end

%% Compute action torques in every instant
if strcmp(UvecFlag,'YES')
    u_vector = zeros(tf/h+1,4);
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
                u_vector(i,1:4) = [u_cte',times(i)];
            end
        case 'VAR'
            ls = length(states);
            for i=1:ls
                u_f= u_function(times(i), u, m, sm);
                u_zoh = zoh_function(times(i),u_f);
                u_vector(i,1:4) = [u_zoh',times(i)];
            end
        case 'CTC_LQR'
            switch goal_Flag
                case 'FIX'
                    ls = length(states);
                    for i=1:ls
                        MFIKsout = sm.MFIKmatrix(states(i,3),states(i,6));
                        M_barsout = sm.M_barmatrix(states(i,3),states(i,6));
                        C_barsout = sm.C_barmatrix(states(i,3),states(i,9),states(i,6),states(i,12));

                        pdiffout(1:3,1) = cp.pss(1:3,1) - states(i,1:3)';
                        pdiffout(4:6,1) = cp.pss(4:6,1) - states(i,7:9)';

                        u_2out = cp.K*pdiffout;

                        uout = (MFIKsout.')*(M_barsout*u_2out + C_barsout*states(i,7:9)');
                        u_zoh = zoh_function(times(i),uout);
                        u_vector(i,1:4) = [u_zoh',times(i)];
                    end
                case 'CIRCLE'
                    ls = length(states);
                    for i=1:ls
                        MFIKsout = sm.MFIKmatrix(states(i,3),states(i,6));
                        M_barsout = sm.M_barmatrix(states(i,3),states(i,6));
                        C_barsout = sm.C_barmatrix(states(i,3),states(i,9),states(i,6),states(i,12));

                        [xssout,yssout,alphassout,xdotssout,ydotssout,alphadotssout,xdotdotssout,ydotdotssout,alphadotdotssout] = TD_circle_of_time_XYA(times(i));
                        pssout = [xssout;yssout;alphassout;xdotssout;ydotssout;alphadotssout];
                        pd_dotdotout = [xdotdotssout; ydotdotssout; alphadotdotssout];

                        pdiffout(1:3,1) = pssout(1:3,1) - states(i,1:3)';
                        pdiffout(4:6,1) = pssout(4:6,1) - states(i,7:9)';

                        u_2out = pd_dotdotout + cp.K*pdiffout;

                        uout = (MFIKsout.')*(M_barsout*u_2out + C_barsout*states(i,7:9)');
                        u_zoh = zoh_function(times(i),uout);
                        u_vector(i,1:4) = [u_zoh',times(i)];
                    end
                case 'POLILINE'
                    ls = length(states);
                    for i=1:ls
                        MFIKsout = sm.MFIKmatrix(states(i,3),states(i,6));
                        M_barsout = sm.M_barmatrix(states(i,3),states(i,6));
                        C_barsout = sm.C_barmatrix(states(i,3),states(i,9),states(i,6),states(i,12));

                        [xssout,yssout,alphassout,xdotssout,ydotssout,alphadotssout,xdotdotssout,ydotdotssout,alphadotdotssout] = TD_polyline_of_time(times(i));
                        pssout = [xssout;yssout;alphassout;xdotssout;ydotssout;alphadotssout];
                        pd_dotdotout = [xdotdotssout; ydotdotssout; alphadotdotssout];

                        pdiffout(1:3,1) = pssout(1:3,1) - states(i,1:3)';
                        pdiffout(4:6,1) = pssout(4:6,1) - states(i,7:9)';

                        u_2out = pd_dotdotout + cp.K*pdiffout;
                        uout = (MFIKsout.')*(M_barsout*u_2out + C_barsout*states(i,7:9)');
                        u_zoh = zoh_function(times(i),uout);
                        u_vector(i,1:4) = [u_zoh',times(i)];
                    end
                otherwise
                    disp('This goal_Flag does not exist for the torques computation')
                    disp('Simulating with a FIX point X=2 Y=2 Alpha = 0')
                    cp.pss = [2,2,0,0,0,0]';
                    ls = length(states);
                    for i=1:ls
                        MFIKsout = sm.MFIKmatrix(states(i,3),states(i,6));
                        M_barsout = sm.M_barmatrix(states(i,3),states(i,6));
                        C_barsout = sm.C_barmatrix(states(i,3),states(i,9),states(i,6),states(i,12));

                        pdiffout(1:3,1) = cp.pss(1:3,1) - states(i,1:3)';
                        pdiffout(4:6,1) = cp.pss(4:6,1) - states(i,7:9)';

                        u_2out = cp.K*pdiffout;

                        uout = (MFIKsout.')*(M_barsout*u_2out + C_barsout*states(i,7:9)');
                        u_zoh = zoh_function(times(i),uout);
                        u_vector(i,1:4) = [u_zoh',times(i)];
                    end
            end 
         case 'CTC_PP'
            switch goal_Flag
                case 'FIX'
                    ls = length(states);
                    for i=1:ls
                        MFIKsout = sm.MFIKmatrix(states(i,3),states(i,6));
                        M_barsout = sm.M_barmatrix(states(i,3),states(i,6));
                        C_barsout = sm.C_barmatrix(states(i,3),states(i,9),states(i,6),states(i,12));

                        pdiffout(1:3,1) = cp.pss(1:3,1) - states(i,1:3)';
                        pdiffout(4:6,1) = cp.pss(4:6,1) - states(i,7:9)';
                        
                        if strcmp(u_IntFlag,'YES')
                            I_kout = Int_Vecout(1,1:3)' + (times(i)-Int_Vecout(1,7))/2*(Int_Vecout(1,4:6)' + pdiffout(1:3,1));
                            u_2out = cp.K*[I_kout;pdiffout];
                            uout = (MFIKsout.')*(M_barsout*u_2out + C_barsout*states(i,7:9)');
                            u_zoh = zoh_function(times(i),uout);
                            u_vector(i,1:4) = [u_zoh',times(i)];
                            
                            Int_Vecout(1,1:3) = I_kout';          % Updating Integral value
                            Int_Vecout(1,4:6) = pdiffout(1:3,1)'; % Updating error value
                            Int_Vecout(1,7) = times(i);           % Updating time value
                        elseif strcmp(u_IntFlag,'NO')
                            u_2out = cp.K*pdiffout;
                            uout = (MFIKsout.')*(M_barsout*u_2out + C_barsout*states(i,7:9)');
                            u_zoh = zoh_function(times(i),uout);
                            u_vector(i,1:4) = [u_zoh',times(i)];
                        else
                            disp('This u_IntFlag does not exist runing the simulation with standart PD controller')
                            u_2out = cp.K*pdiffout;
                            uout = (MFIKsout.')*(M_barsout*u_2out + C_barsout*states(i,7:9)');
                            u_zoh = zoh_function(times(i),uout);
                            u_vector(i,1:4) = [u_zoh',times(i)];
                        end

%                         u_2out = cp.K*pdiffout;
% 
%                         uout = (MFIKsout.')*(M_barsout*u_2out + C_barsout*states(i,7:9)');
%                         u_vector(i,1:4) = [uout',times(i)];
                    end
                case 'CIRCLE'
                    ls = length(states);
                    for i=1:ls
                        MFIKsout = sm.MFIKmatrix(states(i,3),states(i,6));
                        M_barsout = sm.M_barmatrix(states(i,3),states(i,6));
                        C_barsout = sm.C_barmatrix(states(i,3),states(i,9),states(i,6),states(i,12));

                        [xssout,yssout,alphassout,xdotssout,ydotssout,alphadotssout,xdotdotssout,ydotdotssout,alphadotdotssout] = TD_circle_of_time_XYA(times(i));
                        pssout = [xssout;yssout;alphassout;xdotssout;ydotssout;alphadotssout];
                        pd_dotdotout = [xdotdotssout; ydotdotssout; alphadotdotssout];

                        pdiffout(1:3,1) = pssout(1:3,1) - states(i,1:3)';
                        pdiffout(4:6,1) = pssout(4:6,1) - states(i,7:9)';
                        
                        if strcmp(u_IntFlag,'YES')
                            I_kout = Int_Vecout(1,1:3)' + (times(i)-Int_Vecout(1,7))/2*(Int_Vecout(1,4:6)' + pdiffout(1:3,1));
                            u_2out = pd_dotdotout + cp.K*[I_kout;pdiffout];
                            uout = (MFIKsout.')*(M_barsout*u_2out + C_barsout*states(i,7:9)');
                            u_zoh = zoh_function(times(i),uout);
                            u_vector(i,1:4) = [u_zoh',times(i)];
                            
                            Int_Vecout(1,1:3) = I_kout';          % Updating Integral value
                            Int_Vecout(1,4:6) = pdiffout(1:3,1)'; % Updating error value
                            Int_Vecout(1,7) = times(i);           % Updating time value
                        elseif strcmp(u_IntFlag,'NO')
                            u_2out = pd_dotdotout + cp.K*pdiffout;
                            uout = (MFIKsout.')*(M_barsout*u_2out + C_barsout*states(i,7:9)');
                            u_zoh = zoh_function(times(i),uout);
                            u_vector(i,1:4) = [u_zoh',times(i)];
                        else
                            disp('This u_IntFlag does not exist runing the simulation with standart PD controller')
                            u_2out = pd_dotdotout + cp.K*pdiffout;
                            uout = (MFIKsout.')*(M_barsout*u_2out + C_barsout*states(i,7:9)');
                            u_zoh = zoh_function(times(i),uout);
                            u_vector(i,1:4) = [u_zoh',times(i)];
                        end

%                         u_2out = pd_dotdotout + cp.K*pdiffout;
% 
%                         uout = (MFIKsout.')*(M_barsout*u_2out + C_barsout*states(i,7:9)');
%                         u_vector(i,1:4) = [uout',times(i)];
                    end
                case 'POLILINE'
                    ls = length(states);
                    for i=1:ls
                        MFIKsout = sm.MFIKmatrix(states(i,3),states(i,6));
                        M_barsout = sm.M_barmatrix(states(i,3),states(i,6));
                        C_barsout = sm.C_barmatrix(states(i,3),states(i,9),states(i,6),states(i,12));

                        [xssout,yssout,alphassout,xdotssout,ydotssout,alphadotssout,xdotdotssout,ydotdotssout,alphadotdotssout] = TD_polyline_of_time(times(i));
                        pssout = [xssout;yssout;alphassout;xdotssout;ydotssout;alphadotssout];
                        pd_dotdotout = [xdotdotssout; ydotdotssout; alphadotdotssout];

                        pdiffout(1:3,1) = pssout(1:3,1) - states(i,1:3)';
                        pdiffout(4:6,1) = pssout(4:6,1) - states(i,7:9)';
                        
                        if strcmp(u_IntFlag,'YES')
                            I_kout = Int_Vecout(1,1:3)' + (times(i)-Int_Vecout(1,7))/2*(Int_Vecout(1,4:6)' + pdiffout(1:3,1));
                            u_2out = pd_dotdotout + cp.K*[I_kout;pdiffout];
                            uout = (MFIKsout.')*(M_barsout*u_2out + C_barsout*states(i,7:9)');
                            u_zoh = zoh_function(times(i),uout);
                            u_vector(i,1:4) = [u_zoh',times(i)];
                            
                            Int_Vecout(1,1:3) = I_kout';          % Updating Integral value
                            Int_Vecout(1,4:6) = pdiffout(1:3,1)'; % Updating error value
                            Int_Vecout(1,7) = times(i);           % Updating time value
                        elseif strcmp(u_IntFlag,'NO')
                            u_2out = pd_dotdotout + cp.K*pdiffout;
                            uout = (MFIKsout.')*(M_barsout*u_2out + C_barsout*states(i,7:9)');
                            u_zoh = zoh_function(times(i),uout);
                            u_vector(i,1:4) = [u_zoh',times(i)];
                        else
                            disp('This u_IntFlag does not exist runing the simulation with standart PD controller')
                            u_2out = pd_dotdotout + cp.K*pdiffout;
                            uout = (MFIKsout.')*(M_barsout*u_2out + C_barsout*states(i,7:9)');
                            u_zoh = zoh_function(times(i),uout);
                            u_vector(i,1:4) = [u_zoh',times(i)];
                        end

%                         u_2out = pd_dotdotout + cp.K*pdiffout;
%                         uout = (MFIKsout.')*(M_barsout*u_2out + C_barsout*states(i,7:9)');
%                         u_vector(i,1:4) = [uout',times(i)];
                    end
                otherwise
                    disp('This goal_Flag does not exist for the torques computation')
                    disp('Simulating with a FIX point X=2 Y=2 Alpha = 0 with standart PD CTC')
                    cp.pss = [2,2,0,0,0,0]';
                    ls = length(states);
                    for i=1:ls
                        MFIKsout = sm.MFIKmatrix(states(i,3),states(i,6));
                        M_barsout = sm.M_barmatrix(states(i,3),states(i,6));
                        C_barsout = sm.C_barmatrix(states(i,3),states(i,9),states(i,6),states(i,12));

                        pdiffout(1:3,1) = cp.pss(1:3,1) - states(i,1:3)';
                        pdiffout(4:6,1) = cp.pss(4:6,1) - states(i,7:9)';
                        u_2out = cp.K*pdiffout;

                        uout = (MFIKsout.')*(M_barsout*u_2out + C_barsout*states(i,7:9)');
                        u_zoh = zoh_function(times(i),uout);
                        u_vector(i,1:4) = [u_zoh',times(i)];
                    end
            end
        otherwise
            disp('Warning this u_Flag does not exits omiting torques camputation')
            disp('Computation with CTE torques will be launched')
            ls = length(states);
            for i=1:ls
                u_vector(i,1:4) = [u_cte',times(i)];
            end
    end
else
    disp('Plots of the evolution of actions u will not be displayed')
end

%% Compute accelerations
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
figure;
plot(times,states(:,1))
xlabel('t(s)','Interpreter','latex'),ylabel('$x$(meters)','Interpreter','latex'),title('Configuration','Interpreter','latex')

figure;
plot(times,states(:,2))
xlabel('t(s)','Interpreter','latex'),ylabel('$y$(meters)','Interpreter','latex'),title('Configuration','Interpreter','latex')

figure;
plot(times,states(:,3))
xlabel('t(s)','Interpreter','latex'),ylabel('$\alpha$(radians)','Interpreter','latex'),title('Configuration','Interpreter','latex')

figure;
plot(times,states(:,4))
xlabel('t(s)','Interpreter','latex'),ylabel('$\varphi_r$(radians)','Interpreter','latex'),title('Configuration','Interpreter','latex')

figure;
plot(times,states(:,5))
xlabel('t(s)','Interpreter','latex'),ylabel('$\varphi_l$(radians)','Interpreter','latex'),title('Configuration','Interpreter','latex')

figure;
plot(times,states(:,6))
xlabel('t(s)','Interpreter','latex'),ylabel('$\varphi_p$(radians)','Interpreter','latex'),title('Configuration','Interpreter','latex')

figure;
plot(times,states(:,7))
xlabel('t(s)','Interpreter','latex'),ylabel('$\dot{x}$(meters/second)','Interpreter','latex'),title('Velocity','Interpreter','latex')

figure;
plot(times,states(:,8))
xlabel('t(s)','Interpreter','latex'),ylabel('$\dot{y}$(meters/second)','Interpreter','latex'),title('Velocity','Interpreter','latex')

figure;
plot(times,states(:,9))
xlabel('t(s)','Interpreter','latex'),ylabel('$\dot{\alpha}$(radians/second)','Interpreter','latex'),title('Velocity','Interpreter','latex')

figure;
plot(times,states(:,10))
xlabel('t(s)','Interpreter','latex'),ylabel('$\dot{\varphi}_r$(radians/second)','Interpreter','latex'),title('Velocity','Interpreter','latex')

figure;
plot(times,states(:,11))
xlabel('t(s)','Interpreter','latex'),ylabel('$\dot{\varphi}_l$(radians/second)','Interpreter','latex'),title('Velocity','Interpreter','latex')

figure;
plot(times,states(:,12))
xlabel('t(s)','Interpreter','latex'),ylabel('$\dot{\varphi}_p$(radians/second)','Interpreter','latex'),title('Velocity','Interpreter','latex')

figure;
plot(times,qdotdot_out(:,1))
xlabel('t(s)','Interpreter','latex'),ylabel('$\ddot{x}$(meters/$s^2$)','Interpreter','latex'),title('Acceleration','Interpreter','latex')

figure;
plot(times,qdotdot_out(:,2))
xlabel('t(s)','Interpreter','latex'),ylabel('$\ddot{y}$(meters/$s^2$)','Interpreter','latex'),title('Acceleration','Interpreter','latex')

% Plot the action vectors
if strcmp(UvecFlag,'YES')
    figure;
    plot(times,u_vector(:,1:2))
    xlabel('t(s)','Interpreter','latex'),ylabel('$\tau$ of wheels(N$\cdot$m)','Interpreter','latex'),title('Action Torques','Interpreter','latex')
    legend('$\tau$ right','$\tau$ left','Interpreter','latex')

    figure;
    plot(times,u_vector(:,3))
    xlabel('t(s)','Interpreter','latex'),ylabel('$\tau$ pivot(N$\cdot$m)','Interpreter','latex'),title('Action Torques','Interpreter','latex')
end

if strcmp(KinEnFlag,'YES')
    figure;
    plot(times,kinenvalues)
    xlabel('t(s)','Interpreter','latex'),ylabel('Kinetic Energy (Joules)','Interpreter','latex'), title('System Energy','Interpreter','latex')
end

if strcmp(ErrorFlag,'YES') && uswitch
    figure;
    plot(times,errorstates(:,1))
    xlabel('t(s)','Interpreter','latex'),ylabel('$x$ error(meters)','Interpreter','latex'),title('Configuration Error','Interpreter','latex')

    figure;
    plot(times,errorstates(:,2))
    xlabel('t(s)','Interpreter','latex'),ylabel('$y$ error(meters)','Interpreter','latex'),title('Configuration Error','Interpreter','latex')

    figure;
    plot(times,errorstates(:,3))
    xlabel('t(s)','Interpreter','latex'),ylabel('$\alpha$ error(radians)','Interpreter','latex'),title('Configuration Error','Interpreter','latex')

    figure;
    plot(times,errorstates(:,4))
    xlabel('t(s)','Interpreter','latex'),ylabel('$\dot{x}$ error(meters/s)','Interpreter','latex'),title('Velocity Error','Interpreter','latex')

    figure;
    plot(times,errorstates(:,5))
    xlabel('t(s)','Interpreter','latex'),ylabel('$\dot{y}$ error(meters/s)','Interpreter','latex'),title('Velocity Error','Interpreter','latex')

    figure;
    plot(times,errorstates(:,6))
    xlabel('t(s)','Interpreter','latex'),ylabel('$\dot{\alpha}$ error(radians/s)','Interpreter','latex'),title('Velocity Error','Interpreter','latex')
else
    disp('Error plots will not be displayed')
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
                
                TEFlagcheck = strcmp(TargetFlag,'YES') && TEFlag == 1;
                if strcmp(TargetFlag,'NO') || TEFlagcheck
                    for i = 1:n:length(times)
                            q.x = states(i,1);
                            q.y = states(i,2);
                            q.alpha = states(i,3);
                            q.varphi_r = states(i,4);
                            q.varphi_l = states(i,5);
                            q.varphi_p = states(i,6);

                            elapsed = cputime-time;
                            if elapsed > 1200
                                disp(elapsed);
                                disp('took too long to generate the video')
                                break
                            else
                                draw_otbot(m,q) % Drawing frame
                                hold on
                                if strcmp(DistFlag,'YES') && norm(Dist_funct(times(i))) ~= 0
                                    Dist_video = Dist_funct(times(i));
                                    D_XY = Dist_video(1:2,1);
                                    D_XY = D_XY./norm(D_XY);
                                    D_XY = 3*m.l_1*D_XY;
                                    mArrow3([states(i,1) - D_XY(1,1), states(i,2) - D_XY(2,1), 0],[states(i,1), states(i,2), 0], 'color', 'red', 'stemWidth', 3*m.l_1/50);
                                end
                                hold off
                                
                                frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
                                writeVideo(writerObj, frame); 
                            end

                    end
                    close(writerObj); % Saves the movie.
                    
                elseif strcmp(TargetFlag,'YES') && TEFlag == 0
                    for i = 1:n:length(times)
                            q.x = states(i,1);
                            q.y = states(i,2);
                            q.alpha = states(i,3);
                            q.varphi_r = states(i,4);
                            q.varphi_l = states(i,5);
                            q.varphi_p = states(i,6);

                            elapsed = cputime-time;
                            if elapsed > 1200
                                disp(elapsed);
                                disp('took too long to generate the video')
                                break
                            else
                                draw_otbot(m,q)
                                hold on
                                otbot_circle_v2(targetpos(i,1),targetpos(i,2),m.l_2/4,'r');
                                hold off
                                
                                hold on
                                if strcmp(DistFlag,'YES') && norm(Dist_funct(times(i))) ~= 0
                                    Dist_video = Dist_funct(times(i));
                                    D_XY = Dist_video(1:2,1);
                                    D_XY = D_XY./norm(D_XY);
                                    D_XY = 3*m.l_1*D_XY;
                                    mArrow3([states(i,1) - D_XY(1,1), states(i,2) - D_XY(2,1), 0],[states(i,1), states(i,2), 0], 'color', 'red', 'stemWidth', 3*m.l_1/50);
                                end
                                hold off
                                
                                frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
                                writeVideo(writerObj, frame); 
                            end

                    end
                    close(writerObj); % Saves the movie.
                    
                end
                
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
                
                TEFlagcheck = strcmp(TargetFlag,'YES') && TEFlag == 1;
                if strcmp(TargetFlag,'NO') || TEFlagcheck
                    for i = 1:n:length(times)
                            q.x = states(i,1);
                            q.y = states(i,2);
                            q.alpha = states(i,3);
                            q.varphi_r = states(i,4);
                            q.varphi_l = states(i,5);
                            q.varphi_p = states(i,6);

                            elapsed = cputime-time;
                            if elapsed > 1200
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
                                
                                hold on
                                if strcmp(DistFlag,'YES') && norm(Dist_funct(times(i))) ~= 0
                                    Dist_video = Dist_funct(times(i));
                                    D_XY = Dist_video(1:2,1);
                                    D_XY = D_XY./norm(D_XY);
                                    D_XY = 3*m.l_1*D_XY;
                                    mArrow3([states(i,1) - D_XY(1,1), states(i,2) - D_XY(2,1), 0],[states(i,1), states(i,2), 0], 'color', 'red', 'stemWidth', 3*m.l_1/50);
                                end
                                hold off
                                
                                frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
                                writeVideo(writerObj, frame); 
                            end
                    end
                    close(writerObj); % Saves the movie.
                    
                elseif strcmp(TargetFlag,'YES') && TEFlag == 0
                    for i = 1:n:length(times)
                            q.x = states(i,1);
                            q.y = states(i,2);
                            q.alpha = states(i,3);
                            q.varphi_r = states(i,4);
                            q.varphi_l = states(i,5);
                            q.varphi_p = states(i,6);

                            elapsed = cputime-time;
                            if elapsed > 1200
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
                                    hold on
                                    otbot_circle_v2(targetpos(i,1),targetpos(i,2),m.l_2/4,'r');
                                    hold off
                                    jaux2 = jaux2 + 1;
                                    jaux = jaux+1;
                                elseif jaux2>= Gpf
                                    multi_circles(GBpplot,m.l_2/6,'b')
                                    hold on
                                    multi_circles(GPpplot,m.l_2/6,'g')
                                    hold on
                                    otbot_circle_v2(targetpos(i,1),targetpos(i,2),m.l_2/4,'r');
                                    hold off
                                    jaux2 = 1;
                                else
                                    multi_circles(GBpplot,m.l_2/6,'b')
                                    hold on
                                    multi_circles(GPpplot,m.l_2/6,'g')
                                    hold on
                                    otbot_circle_v2(targetpos(i,1),targetpos(i,2),m.l_2/4,'r');
                                    hold off
                                    jaux2 = jaux2 + 1;
                                end
                                
                                hold on
                                if strcmp(DistFlag,'YES') && norm(Dist_funct(times(i))) ~= 0
                                    Dist_video = Dist_funct(times(i));
                                    D_XY = Dist_video(1:2,1);
                                    D_XY = D_XY./norm(D_XY);
                                    D_XY = 3*m.l_1*D_XY;
                                    mArrow3([states(i,1) - D_XY(1,1), states(i,2) - D_XY(2,1), 0],[states(i,1), states(i,2), 0], 'color', 'red', 'stemWidth', 3*m.l_1/50);
                                end
                                hold off
                                
                                frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
                                writeVideo(writerObj, frame); 
                            end
                    end
                    close(writerObj); % Saves the movie.
                    
                end
                                        
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
                
                TEFlagcheck = strcmp(TargetFlag,'YES') && TEFlag == 1;
                if strcmp(TargetFlag,'NO') || TEFlagcheck
                    for i = 1:n:length(times)
                            q.x = states(i,1);
                            q.y = states(i,2);
                            q.alpha = states(i,3);
                            q.varphi_r = states(i,4);
                            q.varphi_l = states(i,5);
                            q.varphi_p = states(i,6);

                            elapsed = cputime-time;
                            if elapsed > 1200
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
                                
                                hold on
                                if strcmp(DistFlag,'YES') && norm(Dist_funct(times(i))) ~= 0
                                    Dist_video = Dist_funct(times(i));
                                    D_XY = Dist_video(1:2,1);
                                    D_XY = D_XY./norm(D_XY);
                                    D_XY = 3*m.l_1*D_XY;
                                    mArrow3([states(i,1) - D_XY(1,1), states(i,2) - D_XY(2,1), 0],[states(i,1), states(i,2), 0], 'color', 'red', 'stemWidth', 3*m.l_1/50);
                                end
                                hold off
                                
                                frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
                                writeVideo(writerObj, frame); 
                            end
                    end
                    close(writerObj); % Saves the movie.
                
                elseif strcmp(TargetFlag,'YES') && TEFlag == 0
                    for i = 1:n:length(times)
                            q.x = states(i,1);
                            q.y = states(i,2);
                            q.alpha = states(i,3);
                            q.varphi_r = states(i,4);
                            q.varphi_l = states(i,5);
                            q.varphi_p = states(i,6);

                            elapsed = cputime-time;
                            if elapsed > 1200
                                disp(elapsed);
                                disp('took too long to generate the video')
                                break
                            else
                                draw_otbot(m,q)
                                if jaux2 == 1
                                    GBpplot = GBpmat(:,1:jaux); 
                                    multi_circles(GBpplot,m.l_2/6,'b')
                                    hold on
                                    otbot_circle_v2(targetpos(i,1),targetpos(i,2),m.l_2/4,'r');
                                    hold off
                                    jaux2 = jaux2 + 1;
                                    jaux = jaux+1;
                                elseif jaux2>= Gpf
                                    multi_circles(GBpplot,m.l_2/6,'b')
                                    hold on
                                    otbot_circle_v2(targetpos(i,1),targetpos(i,2),m.l_2/4,'r');
                                    hold off
                                    jaux2 = 1;
                                else
                                    multi_circles(GBpplot,m.l_2/6,'b')
                                    hold on
                                    otbot_circle_v2(targetpos(i,1),targetpos(i,2),m.l_2/4,'r');
                                    hold off
                                    jaux2 = jaux2 + 1;
                                end
                                
                                hold on
                                if strcmp(DistFlag,'YES') && norm(Dist_funct(times(i))) ~= 0
                                    Dist_video = Dist_funct(times(i));
                                    D_XY = Dist_video(1:2,1);
                                    D_XY = D_XY./norm(D_XY);
                                    D_XY = 3*m.l_1*D_XY;
                                    mArrow3([states(i,1) - D_XY(1,1), states(i,2) - D_XY(2,1), 0],[states(i,1), states(i,2), 0], 'color', 'red', 'stemWidth', 3*m.l_1/50);
                                end
                                hold off
                                
                                frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
                                writeVideo(writerObj, frame); 
                            end
                    end
                    close(writerObj); % Saves the movie.
                end
                
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
                
                TEFlagcheck = strcmp(TargetFlag,'YES') && TEFlag == 1;
                if strcmp(TargetFlag,'NO') || TEFlagcheck
                    for i = 1:n:length(times)
                            q.x = states(i,1);
                            q.y = states(i,2);
                            q.alpha = states(i,3);
                            q.varphi_r = states(i,4);
                            q.varphi_l = states(i,5);
                            q.varphi_p = states(i,6);

                            elapsed = cputime-time;
                            if elapsed > 1200
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
                                
                                hold on
                                if strcmp(DistFlag,'YES') && norm(Dist_funct(times(i))) ~= 0
                                    Dist_video = Dist_funct(times(i));
                                    D_XY = Dist_video(1:2,1);
                                    D_XY = D_XY./norm(D_XY);
                                    D_XY = 3*m.l_1*D_XY;
                                    mArrow3([states(i,1) - D_XY(1,1), states(i,2) - D_XY(2,1), 0],[states(i,1), states(i,2), 0], 'color', 'red', 'stemWidth', 3*m.l_1/50);
                                end
                                hold off
                                
                                frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
                                writeVideo(writerObj, frame); 
                            end
                    end
                    close(writerObj); % Saves the movie.
                
                elseif strcmp(TargetFlag,'YES') && TEFlag == 0
                    for i = 1:n:length(times)
                            q.x = states(i,1);
                            q.y = states(i,2);
                            q.alpha = states(i,3);
                            q.varphi_r = states(i,4);
                            q.varphi_l = states(i,5);
                            q.varphi_p = states(i,6);

                            elapsed = cputime-time;
                            if elapsed > 1200
                                disp(elapsed);
                                disp('took too long to generate the video')
                                break
                            else
                                draw_otbot(m,q)
                                if jaux2 == 1
                                    GPpplot = GPpmat(:,1:jaux);
                                    multi_circles(GPpplot,m.l_2/6,'g')
                                    hold on
                                    otbot_circle_v2(targetpos(i,1),targetpos(i,2),m.l_2/4,'r');
                                    hold off
                                    jaux2 = jaux2 + 1;
                                    jaux = jaux+1;
                                elseif jaux2>= Gpf
                                    multi_circles(GPpplot,m.l_2/6,'g')
                                    hold on
                                    otbot_circle_v2(targetpos(i,1),targetpos(i,2),m.l_2/4,'r');
                                    hold off
                                    jaux2 = 1;
                                else
                                    multi_circles(GPpplot,m.l_2/6,'g')
                                    hold on
                                    otbot_circle_v2(targetpos(i,1),targetpos(i,2),m.l_2/4,'r');
                                    hold off
                                    jaux2 = jaux2 + 1;
                                end
                                
                                hold on
                                if strcmp(DistFlag,'YES') && norm(Dist_funct(times(i))) ~= 0
                                    Dist_video = Dist_funct(times(i));
                                    D_XY = Dist_video(1:2,1);
                                    D_XY = D_XY./norm(D_XY);
                                    D_XY = 3*m.l_1*D_XY;
                                    mArrow3([states(i,1) - D_XY(1,1), states(i,2) - D_XY(2,1), 0],[states(i,1), states(i,2), 0], 'color', 'red', 'stemWidth', 3*m.l_1/50);
                                end
                                hold off
                                
                                frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
                                writeVideo(writerObj, frame); 
                            end
                    end
                    close(writerObj); % Saves the movie.
                end
                      
            otherwise
                disp('This TrackFlag does not exist, omitting the video')
        end            
    case 'NO'
        disp('Simulation without video launched');
    otherwise
        disp('This FlagVideo input does not exist');
end

%% Create output data
Upath = userpath;

ls = length(states);
if h2>h
    hts = uint64(h2/h);
    ls = length(states);
    
    % Build data output y including alpha_dot, x_ddot, y_ddot
    % First we add the proportional Bias aka Scaling factor error
    y_out_data1 = add_bias_IMU(Biasalphadot_Flag, Biasxddot_Flag, Biasyddot_Flag, states, IMUxy, hts, ls, IMU);
    % Then we add the noise
    y_out_data = add_noise_IMU(Noisealphadot_Flag, Noisexddot_Flag, Noiseyddot_Flag, Nopt_Flag1, Nopt_Flag2, Nopt_Flag3, y_out_data1, IMU);
    
    % Build action output 
    u_out_data = u_vector(1:hts:ls,1:3);
    
    % Save the real values of the parameters
    real_pars = give_realpars(m,Biasalphadot_Flag, Biasxddot_Flag, Biasyddot_Flag, IMU);

    % Save the values of the state & accelerations to compare it later
    x_out_data = states(1:hts:ls,:);
    qddot_real_data = IMUxy(1:hts:ls,1:2);
   
    savedir1 = strcat(Upath,'\Model_Identification\GBEOtbot1_V2_IMU_AD(Stp3)\');
    savedir2 = strcat(Upath,'\Model_Identification\GBEOtbot1_V2_IMU_AD(Stp3)\Final_Pars_Sim\Otbot1_CTC_W1_DF_IMU\');
    
    save(strcat(savedir1,'Otbot1_Data.mat'),"y_out_data","u_out_data")
    save(strcat(savedir2,'Otbot1_Data.mat'),"y_out_data","u_out_data")
    save(strcat(savedir2,'Otbot1_DataStates.mat'),"x_out_data")
    save(strcat(savedir2,'Otbot1_DataAccel.mat'),"qddot_real_data")
    save(strcat(savedir1,'initial_states.mat'),"xs0")
    save(strcat(savedir1,'real_pars_data.mat'),"real_pars")
else
    disp('Warning h2 has to be always greater than h omiting output vectors')
end

%% Export the Flags
% Create a struc object to save all the flags and export them

Flags.FlagVideo = FlagVideo;
Flags.Noisealphadot_Flag = Noisealphadot_Flag;   
Flags.Noisexddot_Flag = Noisexddot_Flag;                      
Flags.Noiseyddot_Flag = Noiseyddot_Flag;
Flags.Biasalphadot_Flag = Biasalphadot_Flag;              
Flags.Biasxddot_Flag = Biasxddot_Flag;                     
Flags.Biasyddot_Flag = Biasyddot_Flag;             
Flags.u_Flag = u_Flag;               
Flags.u_IntFlag = u_IntFlag;            
Flags.goal_Flag = goal_Flag;            
Flags.FrictFlag = FrictFlag;            
Flags.DistFlag = DistFlag;           
Flags.TrackFlag = TrackFlag;              
Flags.TargetFlag = TargetFlag;           
Flags.UvecFlag = UvecFlag;                         
Flags.KinEnFlag = KinEnFlag;              
Flags.Jdot_qdotFlag = Jdot_qdotFlag;            
Flags.ErrorFlag = ErrorFlag;
Flags.u = u;
Flags.u_cte = u_cte;
Flags.Fix_point = Fix_point;
Flags.ASR = IMU.ASR;
Flags.XSR = IMU.XSR;
Flags.YSR = IMU.YSR;
Flags.NDgy = IMU.NDgy;
Flags.NDaccX = IMU.NDaccX;
Flags.NDaccY = IMU.NDaccY;
 
save(strcat(savedir2,'Flags_Data.mat'),"Flags")

