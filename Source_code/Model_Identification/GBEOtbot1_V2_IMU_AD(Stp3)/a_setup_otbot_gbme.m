%% 1. Load measured data.

% The measured angular displacement data is loaded and saved as data, 
% an iddata object with a sample time of 0.1 seconds. The set command is 
% used to specify data attributes such as the output name,  output unit, 
% and the start time and units of the time vector.
prompt = 'Please introduce the sample time of the data [0.01s]: ';
TSinput = input(prompt);
if isempty(TSinput)
    TSinput = 0.01;
end

% Set the deviation of the initial guess with respect the real values
prompt = 'Please introduce the deviation of the parameters for the initial guess [50%]: ';
TSinputDev = input(prompt);
if isempty(TSinputDev)
    TSinputDev = 50;
end

% IGD = 1-TSinputDev/100; % This was used in the old version where the
% deviation were applyed as a percentage deviation over the real value
IGAD = TSinputDev/100; % This is the percentage of the maxIG deviation that will be applied to the parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Otbot1_Data.mat'); % Notice that this is not otbot data and needs to be changed
data = iddata(y_out_data, u_out_data, TSinput, 'Name', 'Otbot1_Data'); % Output signal of the system y
                                          % Input signal of the system u
                                          % Sample time Ts
                                          
% Specify input and output names, start time and time units
data.InputName = {'Right wheel torque';      ... % u(1).
                    'Left wheel torque';     ... % u(2).
                    'Platform torque'};          % u(3).
data.InputUnit =  {'Nm';'Nm';'Nm'};
data.OutputName = {'alpha_dot';'x_ddot';'y_ddot'};
data.OutputUnit = {'rad/s';'m/s^2';'m/s^2'};
data.Tstart = 0;
data.TimeUnit = 's';

% Specify intersample behaviour for transformations between discrete time
% and continuous time.
data.InterSample = {'zoh';  % Behaviour for the first input u(1) - Right wheel torque 
                    'zoh';  % Behaviour for the second input u(2) - Left wheel torque 
                    'zoh'}; % Behaviour for the third input u(3) - Platform torque 

%------Plot the Data------%
prompt = 'Do you want to plot the data that has been loaded? Y/N [N]: ';
strinput = input(prompt,'s');
if isempty(strinput)
    strinput = 'N';
end

if strcmp(strinput, 'Y') || strcmp(strinput, 'y')
    plottingsampledata
    prompt = 'Do you want to save the plots? Y/N [N]: '; % Saveing the plots
    strinput = input(prompt,'s');
    if isempty(strinput)
        strinput = 'N';
    end

    if strcmp(strinput, 'Y') || strcmp(strinput, 'y')
        save_sampleplots
    elseif strcmp(strinput, 'N') || strcmp(strinput, 'n')
        disp('The plots will not be saved')
    else
        disp('Unrecognised input, the plots will not be saved')
    end
elseif strcmp(strinput, 'N') || strcmp(strinput, 'n')
    disp('The plots of the loaded data will be omitted')
else
    disp('Unrecognised input, the plots of the loaded data will be omitted')
end
%-------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2. Represent the model using an idnlgrey object
load('real_pars_data.mat')

% Create a nonlinear grey-box model associated with the otbot_m function.
I_b = real_pars(1) + IGmax.Ib*(IGAD*double(~FixorNot(1,1)) + double(FixorNot(1,1)));   % Central moment of inertia of the chassis body about axis 3 [kg*m^2] 
% I_p = real_pars(2) + IGmax.Ip*(IGAD*double(~FixorNot(2,1)) + double(FixorNot(2,1))); % Central moment of inertia of the platform body about axis 3" [kg*m^2]
I_a = real_pars(3) + IGmax.Ia*(IGAD*double(~FixorNot(3,1)) + double(FixorNot(3,1))); % Axial moment of inertia of one wheel [kg*m^2]
I_t = real_pars(4) + IGmax.It*(IGAD*double(~FixorNot(4,1)) + double(FixorNot(4,1))); % Twisting moment of inertia of one wheel [kg*m^2] 

l_1 = real_pars(5) + IGmax.l1*(IGAD*double(~FixorNot(5,1)) + double(FixorNot(5,1))); % Pivot offset relative to the wheels axis [m]
l_2 = real_pars(6) + IGmax.l2*(IGAD*double(~FixorNot(6,1)) + double(FixorNot(6,1))); % One half of the wheels separation [m]

m_b = real_pars(7)  + IGmax.mb*(IGAD*double(~FixorNot(7,1)) + double(FixorNot(7,1))); % Mass of the chassis base [kg]
m_w = real_pars(8)  + IGmax.mw*(IGAD*double(~FixorNot(8,1)) + double(FixorNot(8,1))); % Mass of one wheel [kg] 
m_p = real_pars(9)  + IGmax.mp*(IGAD*double(~FixorNot(9,1)) + double(FixorNot(9,1))); % Mass of the platform [kg]

x_G = real_pars(10) + IGmax.x_G*(IGAD*double(~FixorNot(10,1)) + double(FixorNot(10,1))); % x coord of the c.o.m. of the chassis body in the chassis frame [m]
y_G = real_pars(11) + IGmax.y_G*(IGAD*double(~FixorNot(11,1)) + double(FixorNot(11,1)));     % y coord of the c.o.m. of the chassis body in the chassis frame [m]

x_F = real_pars(12) + IGmax.x_F*(IGAD*double(~FixorNot(12,1)) + double(FixorNot(12,1))); % x coord of the c.o.m. of the platform body in the platform frame [m]
y_F = real_pars(13) + IGmax.x_F*(IGAD*double(~FixorNot(13,1)) + double(FixorNot(13,1))); % y coord of the c.o.m. of the platform body in the platform frame [m]


r = real_pars(14) + IGmax.r*(IGAD*double(~FixorNot(14,1)) + double(FixorNot(14,1)));   % Wheel radius [m]

Abias = real_pars(15) + IGmax.Abias*(IGAD*double(~FixorNot(15,1)) + double(FixorNot(15,1)));   % Bias of alphadot [nounit]

Xbias = real_pars(16) + IGmax.Xbias*(IGAD*double(~FixorNot(16,1)) + double(FixorNot(16,1)));   % Bias of x_ddot [nounit]

Ybias = real_pars(17) + IGmax.Ybias*(IGAD*double(~FixorNot(17,1)) + double(FixorNot(17,1)));   % Bias of y_ddot [nounit]

% Computing the new moment of inertia of the platform as a function of the
% other values x_F,y_F, m_p

I_p = real_pars(2) + m_p*(x_F^2 + y_F^2)*(double(~FixorNot(2,1)) + double(FixorNot(2,1))); % Central moment of inertia of the platform body about axis 3" [kg*m^2]



FileName       = 'otbot1_m_V2_IMU';        % File describing the model structure.
order          = [3 3 12];         % Order — Number of outputs, inputs, and states of the model [Ny,Nu,Nx]
parameters     = {I_b,I_p,I_a,I_t,l_1,l_2,m_b,m_w,m_p,x_G,y_G,x_F,y_F,r,Abias,Xbias,Ybias};  % Parameters of the model

%--------Initial States--------%
prompt = 'HOW DO YOU WANT TO SET THE INITAL STATE? 0/E [E]: ';
mainstrinput = input(prompt,'s');
if isempty(mainstrinput)
    mainstrinput = 'E';
end

if strcmp(mainstrinput, '0') || strcmp(mainstrinput, 'zero')
    initial_states = zeros(12,1); 
elseif strcmp(mainstrinput, 'E') || strcmp(mainstrinput, 'e')
    disp('The initial state will be set to the values of initial_states.mat')
    load('initial_states.mat')
    initial_states = xs0;
else
    disp('Unrecognised input, initial states will be set to 0')
    initial_states = zeros(12,1);
end
%------------------------------%

Ts             = 0;                 % Specify as continuous system.
nlgrm = idnlgrey(FileName,order,parameters,initial_states,Ts,'Name','Otbot1_V2_IMU-GM');

% Specify input and output names, and units.
nlgrm.InputName  = {'Right wheel torque';    ... % u(1).
                    'Left wheel torque';     ... % u(2).
                    'Platform torque'};          % u(3).
                
nlgrm.InputUnit  = {'Nm';'Nm';'Nm'};

nlgrm.OutputName = {'alpha_dot';'x_ddot';'y_ddot'};
                
nlgrm.OutputUnit = {'rad/s';'m/s^2';'m/s^2'};

set(nlgrm,'TimeUnit', 's');

% Specify names and units of the initial states and parameters
% Initial States names and units
nlgrm = setinit(nlgrm, 'Name', {'x';'y';'alpha';'varphi_r';'varphi_l';'varphi_p'; ...
                               'x_dot';'y_dot';'alpha_dot';'varphi_dot_r';'varphi_dot_l';'varphi_dot_p'});

nlgrm = setinit(nlgrm, 'Unit', {'m';'m';'rad';'rad';'rad';'rad'; ...
                                'm/s';'m/s';'rad/s';'rad/s';'rad/s';'rad/s'});

% Inital parameters names and units 
nlgrm = setpar(nlgrm, 'Name', {'Inertia of the chassis body';               ... % I_b
                               'Inertia of the platform body';              ... % I_p
                               'Axial moment of inertia of one wheel';      ... % I_a
                               'Twisting moment of inertia of one wheel';   ... % I_t
                               'Pivot offset';                              ... % l_1
                               'One half of the wheels separation';         ... % l_2
                               'Mass of the chassis base';                  ... % m_b
                               'Mass of one wheel';                         ... % m_w
                               'Mass of the platform';                      ... % m_p
                               'x coord of the c.o.m. of the chassis body'; ... % x_G
                               'y coord of the c.o.m. of the chassis body'; ... % y_G
                               'x coord of the c.o.m. of the plat. body';   ... % x_F
                               'y coord of the c.o.m. of the plat. body';   ... % y_F
                               'Wheel radius';                              ... % r
                               'Alphadot bias';                             ... % Abias 
                               'x_ddot bias';                               ... % Xbias
                               'y_ddot bias'});                                 % Ybias                                                             
                               
nlgrm = setpar(nlgrm, 'Unit', {'kg*m^2';'kg*m^2';'kg*m^2';'kg*m^2'; ...
                                'm';'m';'kg';'kg';'kg';'m';'m';'m';'m';'m';'nounit';'nounit';'nounit'});
                                
nlgrm = setpar(nlgrm, 'Minimum', num2cell(  [eps(0)*ones(4,1);  ... % Moments of Inertia are > 0!
                                             eps(0)*ones(2,1);  ... % Lenghts are > 0 
                                             eps(0)*ones(3,1);  ... % Masses are  > 0
                                             -100*ones(4,1);    ... % We set a min to the coordinates to use fmincon
                                             eps(0);            ... % The value of the wheel radius is > 0!
                                             -100;              ... % Minimum value for Abias not very rellevant
                                             -100;              ... % Minimum value for Xbias not very rellevant
                                             -100]));               % Minimum value for Ybias not very rellevant
                                         
nlgrm = setpar(nlgrm, 'Maximum', num2cell(  [200*ones(4,1);  ... % Moments of Inertia arbitrari maximum to use fmincon
                                             100*ones(2,1);  ... % Lenghts arbitrari maximum to use fmincon
                                             1e6*ones(3,1);  ... % Masses arbitrari maximum to use fmincon
                                             100*ones(4,1);  ... % We set a maximum to the coordinates to use fmincon
                                             30;             ... % Set a maximum value for the wheel radius to use fmincon
                                             100;            ... % We now set a maximum for Abias even that it won't be very important
                                             100;            ... % We now set a maximum for Xbias even that it won't be very important
                                             100]));             % We now set a maximum for Ybias even that it won't be very important             


