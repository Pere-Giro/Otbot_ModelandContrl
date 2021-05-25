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

IGD = 1-TSinputDev/100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('MotorShaft_Data.mat'); % Notice that this is not otbot data and needs to be changed
data = iddata(y_out_data, u_out_data, TSinput, 'Name', 'MotorShaft_Data'); % Output signal of the system y
                                          % Input signal of the system u
                                          % Sample time Ts
                                          
% Specify input and output names, start time and time units
data.InputName = {'torque'};          % u.
data.InputUnit =  {'Nm'};
data.OutputName = {'varphi_dot'};
data.OutputUnit = {'rad/s'};
data.Tstart = 0;
data.TimeUnit = 's';

% Specify intersample behaviour for transformations between discrete time
% and continuous time.
data.InterSample = {'zoh'}; % Behaviour for the input u - MotorTorque 

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
I = real_pars(1) *(IGD*double(~FixorNot(1,1)) + double(FixorNot(1,1))); % Moment of inertia of the motorshaft [kg*m^2]
b = real_pars(2) *(IGD*double(~FixorNot(2,1)) + double(FixorNot(2,1))); % Viscous friction coefficient [kg*m^2*s^-1]


FileName       = 'motorshaft';  % File describing the model structure.
order          = [1 1 2];       % Order — Number of outputs, inputs, and states of the model [Ny,Nu,Nx]
parameters     = {I,b};         % Parameters of the model


%--------Initial States--------%
prompt = 'HOW DO YOU WANT TO SET THE INITAL STATE? 0/E [E]: ';
mainstrinput = input(prompt,'s');
if isempty(mainstrinput)
    mainstrinput = 'E';
end

if strcmp(mainstrinput, '0') || strcmp(mainstrinput, 'zero')
    initial_states = zeros(2,1); 
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
nlgrm = idnlgrey(FileName,order,parameters,initial_states,Ts,'Name','MotorShaft_GM');

% Specify input and output names, and units.
nlgrm.InputName  = {'torque'};          % u
                
nlgrm.InputUnit  = {'Nm'};

nlgrm.OutputName = {'varphi_dot'};
                
nlgrm.OutputUnit = {'rad/s'};

set(nlgrm,'TimeUnit', 's');

% Specify names and units of the initial states and parameters
% Initial States names and units
nlgrm = setinit(nlgrm, 'Name', {'varphi';'varphi_dot'});

nlgrm = setinit(nlgrm, 'Unit', {'rad';'rad/s'});

% Inital parameters names and units 
nlgrm = setpar(nlgrm, 'Name', {'Inertia of the motor shaft';               ... % I
                               'Viscous friction coefficient'});               % b                                                             
                               
nlgrm = setpar(nlgrm, 'Unit', {'kg*m^2';'kg*m^2*s^-1'});
                                
nlgrm = setpar(nlgrm, 'Minimum', num2cell(  [eps(0);  ... % Moments of Inertia are > 0!
                                             -100]));     % Minimum value for viscous friction coefficient
                                         
nlgrm = setpar(nlgrm, 'Maximum', num2cell(  [100;  ... % Moments of Inertia arbitrari maximum to use fmincon
                                             100]));   % We now set a maximum for viscous friction coefficient             


