clearvars
close all
clc
%% Flags
FlagVideo = 'NO'; % Flag to define if you want video output or not
               % (FlagVideo == YES)> simulation with VIDEO output
               % (FlagVideo == NO )> simulation without VIDEO output
               
u_Flag = 'CTC_PP'; % Flag to define if you want a constnat torque for all simulation or not
               % (u_Flag == CTE)> simulation with constant torques
               % (u_Flag == VAR)> simulation with a torque function of time
               % (u_Flag == CTC_LQR)> simulation with a CTC controler using
               % LQR
               % (u_Flag == CTC_PP)> simulation with a CTC controler using
               % pole placement method
                
                % Now we create torques comand vector (Only used if u_Flag is set to CTE or VAR)
                u = [0.1;0.1;-0.1]; % All in Nm
                
u_IntFlag = 'NO'; % Flag to define if you want a controler including an Integrative component or not
               % (u_IntFlag == YES) the CTC controller tuned via PP will
               % include an integrative term.
               % (u_IntFlag == NO) the CTC controller tuned via PP will not
               % include an integrative term.
               
goal_Flag = 'POLILINE'; % Flag to define the desired trajectory for the robot, only works if u is set to LQR or PP
               % (goal_Flag == FIX)> goal trajectory is a fixed point
               % (goal_Flag == CIRCLE)> goal trajectory is a circle
               % (goal_Flag == POLILINE)> goal trajectory is a poliline
               
               % Set a fixed point in case goal_Flag = FIX
                Fix_point = [2,2,0,0,0,0]'; % [x,y,alpha,varphi_r,varphi_l,varphi_p] velocities are 0
               
FrictFlag = 'NO'; % Flag to define if we want to add friction in our model or not
               % (FrictFlag == YES)> Our model will have frictions
               % (FrictFlag == NO)> Our model will not contain frictions
               
DistFlag = 'NO'; % Flag to define if we want to add non modelled disturbances during the simulation
               % (DistFlag == YES)> Non modelled disturbances will be added 
               % (DistFlag == NO)> Non modelled disturbances will NOT be added
               
TrackFlag = 'NO'; % Flag to define if you want to track the path of each center of mass
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
                
TracePath = 'NO'; % Flag to decide how the tracking or the traces of the point will be done
                    % (TracePath == YES)> The traces of the paths will be
                    % done with * and continuous lines
                    % (Tracepath == NO)> The traces of the paths will be
                    % dots or circles with a fixed timestep all along the
                    % trajectory (old way of doing it)
                    colorTP = 'b';
               
UvecFlag = 'YES'; % Flag to choose if you want to compute and display the plots of action vector u
                % (UvecFlag == YES)> action vector u plots will be displayed
                % (UvecFlag == NO)> action vector u plots will not be displayed
                           
KinEnFlag = 'NO'; % Flag to choose if Kinetic energic plot must be displayed or not
               % (KinEnFlag == YES)> Outputs a plot of kinetic energy
               % (KinEnFlag == NO )> No output plot of kinetic energy
               
Jdot_qdotFlag = 'NO'; % Flag to choose if plot of equation Jdot*qdot is equal to 0 during all simulation
               % (Jdot_qdotFlag == YES)> Outputs a plot of this equation
               % (Jdot_qdotFlag == NO)> No output of this plot 
               
ErrorFlag = 'YES'; % Flag to choose if error plots are displayed or not
               % (ErrorFlag == YES)> Error plots will be displayed
               % (ErrorFlag == NO)> Error plots will not be displayed
               
%% Setting up simulation parameters

tf = 30;      % Final simulation time
h = 0.001;    % Number of samples within final vectors (times and states)
opts = odeset('MaxStep',1e-3); % Max integration time step

global tcontrol
tcontrol = [];

%% Setting up model parameters & matrices
load("m_struc")
load("sm_struc")

%% Initial conditions
p0 = [0,0,0]';
varphi0 = [0,0,0]';

% pdot0 will be set in order to fullfill the speeds of the desired
% trajectory in order to have 0 errors in configuration and speeds at time
% 0. This is to perform NoKpKv test.

% pdot0 = [0.3*cos(0.3*0),-0.3*sin(0.3*0),0]';

pdot0 = zeros(3,1);

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
% cp struct object stands for control parameters in contains the desires
% goal fix point and the gain matrix

%% Creting gain matrices for the CTC 

% Our linear system matrices are
Actc = [zeros(3,3),eye(3,3); zeros(3,3), zeros(3,3)];
Bctc = [zeros(3,3);eye(3,3)];

% Defining the poles
s1p = -4/3;
s2p = 10*s1p;

if strcmp(u_Flag,'CTC_LQR')
    Qctc = [3e5*eye(3,3),zeros(3,3);
            zeros(3,3), 5*eye(3,3)];
    Rctc = 5*eye(3,3);
    cp.K = lqr(Actc,Bctc,Qctc,Rctc);
    
elseif strcmp(u_Flag,'CTC_PP')
    if strcmp(u_IntFlag,'NO')
        % EIG_set = [-0.8,-0.8,-0.8, -8, -8, -8];
        EIG_set = [s1p, s1p, s1p, s2p, s2p, s2p];
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
        % EIG_set = [-0.8,-0.8,-0.8, -8, -8, -8];
        EIG_set = [s1p, s1p, s1p, s2p, s2p, s2p];
        cp.K = place(Actc,Bctc,EIG_set);
    end
else
    cp.K = 0;
end

%% Setting-up the Int_Vec (Integrative Vector)
% This part is designed to create a vector that will contain the previous
% values of the Integral and it's time intant. We need this, beacuse we
% will be approximating the integral using the trapezoidal rule.

Int_Vec = zeros(1,7); %[I_k-1(1,3), error_k-1(1,3), t_k-1(1)]

%% Building the diferential equation
dxdt= @(t,xs) xdot2_otbot1_CTC_FDPFD(t, xs, m, sm, u, u_Flag, cp, goal_Flag, FrictFlag, DistFlag, u_IntFlag, Int_Vec);

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
    disp('The combination of flags used is not correct to compute Target')
    disp('Please change the flags to a possible combination in order to display the desired results')
    disp('Target loaction will not be computed')
end

%% Compute action torques in every instant
if strcmp(UvecFlag,'YES')
    u_vector = zeros(tf/h+1,4);
    Int_Vecout = zeros(1,7); %[I_k-1(1,3), error_k-1(1,3), t_k-1(1)]
    switch u_Flag
        case 'CTE'
            ls = length(states);
            for i=1:ls
                u_vector(i,1:4) = [u',times(i)];
            end
        case 'VAR'
            ls = length(states);
            for i=1:ls
                u_f= u_function(times(i),u);
                u_vector(i,1:4) = [u_f',times(i)];
            end
        case 'CTC_LQR'
            switch goal_Flag
                case 'FIX'
                    ls = length(states);
                    for i=1:ls
                        M_barsout2 = sm.M_barmatrix2(states(i,3),states(i,6));
                        C_barsout2 = sm.C_barmatrix2(states(i,3),states(i,9),states(i,6),states(i,12));

                        pdiffout(1:3,1) = cp.pss(1:3,1) - states(i,1:3)';
                        pdiffout(4:6,1) = cp.pss(4:6,1) - states(i,7:9)';

                        u_2out = cp.K*pdiffout;

                        uout = M_barsout2*u_2out + C_barsout2*states(i,7:9)';
                        u_vector(i,1:4) = [uout',times(i)];
                    end
                case 'CIRCLE'
                    ls = length(states);
                    for i=1:ls
                        M_barsout2 = sm.M_barmatrix2(states(i,3),states(i,6));
                        C_barsout2 = sm.C_barmatrix2(states(i,3),states(i,9),states(i,6),states(i,12));

                        [xssout,yssout,alphassout,xdotssout,ydotssout,alphadotssout,xdotdotssout,ydotdotssout,alphadotdotssout] = TD_circle_of_time_XYA(times(i));
                        pssout = [xssout;yssout;alphassout;xdotssout;ydotssout;alphadotssout];
                        pd_dotdotout = [xdotdotssout; ydotdotssout; alphadotdotssout];

                        pdiffout(1:3,1) = pssout(1:3,1) - states(i,1:3)';
                        pdiffout(4:6,1) = pssout(4:6,1) - states(i,7:9)';

                        u_2out = pd_dotdotout + cp.K*pdiffout;

                        uout = M_barsout2*u_2out + C_barsout2*states(i,7:9)';
                        u_vector(i,1:4) = [uout',times(i)];
                    end
                case 'POLILINE'
                    ls = length(states);
                    for i=1:ls
                        M_barsout2 = sm.M_barmatrix2(states(i,3),states(i,6));
                        C_barsout2 = sm.C_barmatrix2(states(i,3),states(i,9),states(i,6),states(i,12));

                        [xssout,yssout,alphassout,xdotssout,ydotssout,alphadotssout,xdotdotssout,ydotdotssout,alphadotdotssout] = TD_polyline_of_time(times(i));
                        pssout = [xssout;yssout;alphassout;xdotssout;ydotssout;alphadotssout];
                        pd_dotdotout = [xdotdotssout; ydotdotssout; alphadotdotssout];

                        pdiffout(1:3,1) = pssout(1:3,1) - states(i,1:3)';
                        pdiffout(4:6,1) = pssout(4:6,1) - states(i,7:9)';

                        u_2out = pd_dotdotout + cp.K*pdiffout;
                        uout = M_barsout2*u_2out + C_barsout2*states(i,7:9)';
                        u_vector(i,1:4) = [uout',times(i)];
                    end
                otherwise
                    disp('This goal_Flag does not exist for the torques computation')
                    disp('Simulating with a FIX point X=2 Y=2 Alpha = 0')
                    cp.pss = [2,2,0,0,0,0]';
                    ls = length(states);
                    for i=1:ls
                        M_barsout2 = sm.M_barmatrix2(states(i,3),states(i,6));
                        C_barsout2 = sm.C_barmatrix2(states(i,3),states(i,9),states(i,6),states(i,12));

                        pdiffout(1:3,1) = cp.pss(1:3,1) - states(i,1:3)';
                        pdiffout(4:6,1) = cp.pss(4:6,1) - states(i,7:9)';

                        u_2out = cp.K*pdiffout;

                        uout = M_barsout2*u_2out + C_barsout2*states(i,7:9)';
                        u_vector(i,1:4) = [uout',times(i)];
                    end
            end 
         case 'CTC_PP'
            switch goal_Flag
                case 'FIX'
                    ls = length(states);
                    for i=1:ls
                        M_barsout2 = sm.M_barmatrix2(states(i,3),states(i,6));
                        C_barsout2 = sm.C_barmatrix2(states(i,3),states(i,9),states(i,6),states(i,12));

                        pdiffout(1:3,1) = cp.pss(1:3,1) - states(i,1:3)';
                        pdiffout(4:6,1) = cp.pss(4:6,1) - states(i,7:9)';
                        
                        if strcmp(u_IntFlag,'YES')
                            I_kout = Int_Vecout(1,1:3)' + (times(i)-Int_Vecout(1,7))/2*(Int_Vecout(1,4:6)' + pdiffout(1:3,1));
                            u_2out = cp.K*[I_kout;pdiffout];
                            uout = M_barsout2*u_2out + C_barsout2*states(i,7:9)';
                            u_vector(i,1:4) = [uout',times(i)];
                            
                            Int_Vecout(1,1:3) = I_kout';          % Updating Integral value
                            Int_Vecout(1,4:6) = pdiffout(1:3,1)'; % Updating error value
                            Int_Vecout(1,7) = times(i);           % Updating time value
                        elseif strcmp(u_IntFlag,'NO')
                            u_2out = cp.K*pdiffout;
                            uout = M_barsout2*u_2out + C_barsout2*states(i,7:9)';        
                            u_vector(i,1:4) = [uout',times(i)];
                        else
                            disp('This u_IntFlag does not exist runing the simulation with standart PD controller')
                            u_2out = cp.K*pdiffout;
                            uout = M_barsout2*u_2out + C_barsout2*states(i,7:9)';        
                            u_vector(i,1:4) = [uout',times(i)];
                        end

%                         u_2out = cp.K*pdiffout;
% 
%                         uout = (MFIKsout.')*(M_barsout*u_2out + C_barsout*states(i,7:9)');
%                         u_vector(i,1:4) = [uout',times(i)];
                    end
                case 'CIRCLE'
                    ls = length(states);
                    for i=1:ls
                        M_barsout2 = sm.M_barmatrix2(states(i,3),states(i,6));
                        C_barsout2 = sm.C_barmatrix2(states(i,3),states(i,9),states(i,6),states(i,12));

                        [xssout,yssout,alphassout,xdotssout,ydotssout,alphadotssout,xdotdotssout,ydotdotssout,alphadotdotssout] = TD_circle_of_time_XYA(times(i));
                        pssout = [xssout;yssout;alphassout;xdotssout;ydotssout;alphadotssout];
                        pd_dotdotout = [xdotdotssout; ydotdotssout; alphadotdotssout];

                        pdiffout(1:3,1) = pssout(1:3,1) - states(i,1:3)';
                        pdiffout(4:6,1) = pssout(4:6,1) - states(i,7:9)';
                        
                        if strcmp(u_IntFlag,'YES')
                            I_kout = Int_Vecout(1,1:3)' + (times(i)-Int_Vecout(1,7))/2*(Int_Vecout(1,4:6)' + pdiffout(1:3,1));
                            u_2out = pd_dotdotout + cp.K*[I_kout;pdiffout];
                            uout = M_barsout2*u_2out + C_barsout2*states(i,7:9)';
                            u_vector(i,1:4) = [uout',times(i)];
                            
                            Int_Vecout(1,1:3) = I_kout';          % Updating Integral value
                            Int_Vecout(1,4:6) = pdiffout(1:3,1)'; % Updating error value
                            Int_Vecout(1,7) = times(i);           % Updating time value
                        elseif strcmp(u_IntFlag,'NO')
                            u_2out = pd_dotdotout + cp.K*pdiffout;
                            uout = M_barsout2*u_2out + C_barsout2*states(i,7:9)';        
                            u_vector(i,1:4) = [uout',times(i)];
                        else
                            disp('This u_IntFlag does not exist runing the simulation with standart PD controller')
                            u_2out = pd_dotdotout + cp.K*pdiffout;
                            uout = M_barsout2*u_2out + C_barsout2*states(i,7:9)';        
                            u_vector(i,1:4) = [uout',times(i)];
                        end

%                         u_2out = pd_dotdotout + cp.K*pdiffout;
% 
%                         uout = (MFIKsout.')*(M_barsout*u_2out + C_barsout*states(i,7:9)');
%                         u_vector(i,1:4) = [uout',times(i)];
                    end
                case 'POLILINE'
                    ls = length(states);
                    for i=1:ls
                        M_barsout2 = sm.M_barmatrix2(states(i,3),states(i,6));
                        C_barsout2 = sm.C_barmatrix2(states(i,3),states(i,9),states(i,6),states(i,12));

                        [xssout,yssout,alphassout,xdotssout,ydotssout,alphadotssout,xdotdotssout,ydotdotssout,alphadotdotssout] = TD_polyline_of_time(times(i));
                        pssout = [xssout;yssout;alphassout;xdotssout;ydotssout;alphadotssout];
                        pd_dotdotout = [xdotdotssout; ydotdotssout; alphadotdotssout];

                        pdiffout(1:3,1) = pssout(1:3,1) - states(i,1:3)';
                        pdiffout(4:6,1) = pssout(4:6,1) - states(i,7:9)';
                        
                        if strcmp(u_IntFlag,'YES')
                            I_kout = Int_Vecout(1,1:3)' + (times(i)-Int_Vecout(1,7))/2*(Int_Vecout(1,4:6)' + pdiffout(1:3,1));
                            u_2out = pd_dotdotout + cp.K*[I_kout;pdiffout];
                            uout = M_barsout2*u_2out + C_barsout2*states(i,7:9)';
                            u_vector(i,1:4) = [uout',times(i)];
                            
                            Int_Vecout(1,1:3) = I_kout';          % Updating Integral value
                            Int_Vecout(1,4:6) = pdiffout(1:3,1)'; % Updating error value
                            Int_Vecout(1,7) = times(i);           % Updating time value
                        elseif strcmp(u_IntFlag,'NO')
                            u_2out = pd_dotdotout + cp.K*pdiffout;
                            uout = M_barsout2*u_2out + C_barsout2*states(i,7:9)';        
                            u_vector(i,1:4) = [uout',times(i)];
                        else
                            disp('This u_IntFlag does not exist runing the simulation with standart PD controller')
                            u_2out = pd_dotdotout + cp.K*pdiffout;
                            uout = M_barsout2*u_2out + C_barsout2*states(i,7:9)';        
                            u_vector(i,1:4) = [uout',times(i)];
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
                        M_barsout2 = sm.M_barmatrix2(states(i,3),states(i,6));
                        C_barsout2 = sm.C_barmatrix2(states(i,3),states(i,9),states(i,6),states(i,12));

                        pdiffout(1:3,1) = cp.pss(1:3,1) - states(i,1:3)';
                        pdiffout(4:6,1) = cp.pss(4:6,1) - states(i,7:9)';
                        u_2out = cp.K*pdiffout;

                        uout = M_barsout2*u_2out + C_barsout2*states(i,7:9)';
                        u_vector(i,1:4) = [uout',times(i)];
                    end
            end
        otherwise
            disp('Warning this u_Flag does not exits omiting torques camputation')
            disp('Computation with CTE torques will be launched')
            ls = length(states);
            for i=1:ls
                u_vector(i,1:4) = [u',times(i)];
            end
    end
else
    disp('Plots of the evolution of actions u will not be displayed')
end

%% Plot results
figure('WindowState','maximized');
plot(times,states(:,1))
xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$x$ [m]','Interpreter','latex','FontSize', 17),title('Configuration','Interpreter','latex','FontSize', 17)
grid on

figure('WindowState','maximized');
plot(times,states(:,2))
xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$y$ [m]','Interpreter','latex','FontSize', 17),title('Configuration','Interpreter','latex','FontSize', 17)
grid on

figure('WindowState','maximized');
plot(times,states(:,3))
xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\alpha$ [rad]','Interpreter','latex','FontSize', 17),title('Configuration','Interpreter','latex','FontSize', 17)
grid on

figure('WindowState','maximized');
plot(times,states(:,4))
xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\varphi_r$ [rad]','Interpreter','latex','FontSize', 17),title('Configuration','Interpreter','latex','FontSize', 17)
grid on

figure('WindowState','maximized');
plot(times,states(:,5))
xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\varphi_l$ [rad]','Interpreter','latex','FontSize', 17),title('Configuration','Interpreter','latex','FontSize', 17)
grid on

figure('WindowState','maximized');
plot(times,states(:,6))
xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\varphi_p$ [rad]','Interpreter','latex','FontSize', 17),title('Configuration','Interpreter','latex','FontSize', 17)
grid on

figure('WindowState','maximized');
plot(times,states(:,7))
xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\dot{x}$ [m/s]','Interpreter','latex','FontSize', 17),title('Velocity','Interpreter','latex','FontSize', 17)
grid on

figure('WindowState','maximized');
plot(times,states(:,8))
xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\dot{y}$ [m/s]','Interpreter','latex','FontSize', 17),title('Velocity','Interpreter','latex','FontSize', 17)
grid on

figure('WindowState','maximized');
plot(times,states(:,9))
xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\dot{\alpha}$ [rad/s]','Interpreter','latex','FontSize', 17),title('Velocity','Interpreter','latex','FontSize', 17)
grid on

figure('WindowState','maximized');
plot(times,states(:,10))
xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\dot{\varphi}_r$ [rad/s]','Interpreter','latex','FontSize', 17),title('Velocity','Interpreter','latex','FontSize', 17)
grid on

figure('WindowState','maximized');
plot(times,states(:,11))
xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\dot{\varphi}_l$ [rad/s]','Interpreter','latex','FontSize', 17),title('Velocity','Interpreter','latex','FontSize', 17)
grid on

figure('WindowState','maximized');
plot(times,states(:,12))
xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\dot{\varphi}_p$ [rad/s]','Interpreter','latex','FontSize', 17),title('Velocity','Interpreter','latex','FontSize', 17)
grid on

% Plot the action vectors
if strcmp(UvecFlag,'YES')
    figure('WindowState','maximized');
    plot(times,u_vector(:,1:2))
    xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\tau$ of wheels [N$\cdot$m]','Interpreter','latex','FontSize', 17),title('Action Torques','Interpreter','latex','FontSize', 17)
    legend('$\tau$ right','$\tau$ left','Interpreter','latex','FontSize', 17)
    grid on
    
    figure('WindowState','maximized');
    plot(times,u_vector(:,3))
    xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\tau$ pivot [N$\cdot$m]','Interpreter','latex','FontSize', 17),title('Action Torques','Interpreter','latex','FontSize', 17)
    grid on
end

if strcmp(KinEnFlag,'YES')
    figure('WindowState','maximized');
    plot(times,kinenvalues)
    xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('Kinetic Energy [J]','Interpreter','latex','FontSize', 17), title('System Energy','Interpreter','latex','FontSize', 17)
    grid on
end

% The errors have to be changed in sign beacause of the new criterion 
if strcmp(ErrorFlag,'YES') && uswitch
    figure('WindowState','maximized');
    plot(times,-errorstates(:,1))
    xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$x$ error [m]','Interpreter','latex','FontSize', 17),title('Configuration Error','Interpreter','latex','FontSize', 17)
    grid on
    
    figure('WindowState','maximized');
    plot(times,-errorstates(:,2))
    xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$y$ error [m]','Interpreter','latex','FontSize', 17),title('Configuration Error','Interpreter','latex','FontSize', 17)
    grid on
    
    figure('WindowState','maximized');
    plot(times,-errorstates(:,3))
    xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\alpha$ error [rad]','Interpreter','latex','FontSize', 17),title('Configuration error in $\alpha$','Interpreter','latex','FontSize', 17)
    grid on
    
    figure('WindowState','maximized');
    plot(times,-errorstates(:,4))
    xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\dot{x}$ error [m/s]','Interpreter','latex','FontSize', 17),title('Velocity Error','Interpreter','latex','FontSize', 17)
    grid on
    
    figure('WindowState','maximized');
    plot(times,-errorstates(:,5))
    xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\dot{y}$ error [m/s]','Interpreter','latex','FontSize', 17),title('Velocity Error','Interpreter','latex','FontSize', 17)
    grid on
    
    figure('WindowState','maximized');
    plot(times,-errorstates(:,6))
    xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\dot{\alpha}$ error [m/s]','Interpreter','latex','FontSize', 17),title('Velocity error in $\dot{\alpha}$','Interpreter','latex','FontSize', 17)
    grid on
    
    
    % Plot error x and y in same plot
    figure('WindowState','maximized');
    plot(times,-errorstates(:,1))
    hold on
    plot(times,-errorstates(:,2))
    hold off
    xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('error [m]','Interpreter','latex','FontSize', 17),title('Configuration error in $x$ and $y$','Interpreter','latex','FontSize', 17)
    legend('$x$ error','$y$ error','Interpreter','latex','FontSize', 17)
    grid on
    
    % Plot error x_dot and y_dot in same plot
    figure('WindowState','maximized');
    plot(times,-errorstates(:,4))
    hold on
    plot(times,-errorstates(:,5))
    hold off
    xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('error [m/s]','Interpreter','latex','FontSize', 17),title('Velocity error in $\dot{x}$ and $\dot{y}$','Interpreter','latex','FontSize', 17)
    legend('$\dot{x}$ error','$\dot{y}$ error','Interpreter','latex','FontSize', 17)
    grid on
    
else
    disp('Error plots will not be displayed')
end

%% Seting Jdot_qdot = 0 Flag

if strcmp(Jdot_qdotFlag,'YES')
    ls = length(states);
    J_qvalues = zeros(ls,3);
    for i=1:ls
        Jplot = sm.Jmatrix(states(i,3),states(i,6));
        qdotplot = states(i,7:12)';
        J_qvalues(i,:) = (Jplot*qdotplot)';
    end
    figure;
    plot(times,J_qvalues)
    xlabel('t(s)'),ylabel('Result vector'), title('Vector result of J*qdot')
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
                                if strcmp(TracePath,'YES')
                                    GPpplotcont = states(1:i,1:2)';
                                    tracepathscont(GPpplotcont,colorTP);
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
                                if strcmp(TracePath,'YES')
                                    GPpplotcont = states(1:i,1:2)';
                                    tracepathscont(GPpplotcont,colorTP);
                                end
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
