%% Specify Initial Sate estimation
% By default the inital states will not be estimated if we want to estimate
% them we must set it up.

nlgrm.InitialStates(1).Fixed = true;  % x          [m]
nlgrm.InitialStates(2).Fixed = true;  % y          [m]
nlgrm.InitialStates(3).Fixed = true;  % alpha      [rad]
nlgrm.InitialStates(4).Fixed = true;  % varphi_r   [rad]
nlgrm.InitialStates(5).Fixed = true;  % varphi_l   [rad]
nlgrm.InitialStates(6).Fixed = true;  % varphi_p   [rad]
nlgrm.InitialStates(7).Fixed = true;  % x_dot      [m/s]
nlgrm.InitialStates(8).Fixed = true;  % y_dot      [m/s]
nlgrm.InitialStates(9).Fixed = true;  % alpha_dot  [rad/s]
nlgrm.InitialStates(10).Fixed = true; % varphi_r   [rad/s]
nlgrm.InitialStates(11).Fixed = true; % varphi_l   [rad/s]
nlgrm.InitialStates(12).Fixed = true; % varphi_p   [rad/s]

%% Specify known parameters
% Specify which parameters are alredy known and which wants we need to
% estimate

nlgrm.parameters(1).Fixed = FixorNot(1,1);  % I_b [kg*m^2]
nlgrm.parameters(2).Fixed = FixorNot(2,1);  % I_p [kg*m^2]
nlgrm.parameters(3).Fixed = FixorNot(3,1);  % I_a [kg*m^2]
nlgrm.parameters(4).Fixed = FixorNot(4,1);  % I_t [kg*m^2]
nlgrm.parameters(5).Fixed = FixorNot(5,1);  % l_1 [m]
nlgrm.parameters(6).Fixed = FixorNot(6,1);  % l_2 [m]
nlgrm.parameters(7).Fixed = FixorNot(7,1);  % m_b [kg] 
nlgrm.parameters(8).Fixed = FixorNot(8,1);  % m_w [kg]
nlgrm.parameters(9).Fixed = FixorNot(9,1);  % m_p [kg] 
nlgrm.parameters(10).Fixed = FixorNot(10,1); % x_G [m]
nlgrm.parameters(11).Fixed = FixorNot(11,1); % y_G [m]
nlgrm.parameters(12).Fixed = FixorNot(12,1); % x_F [m]
nlgrm.parameters(13).Fixed = FixorNot(13,1); % y_F [m]
nlgrm.parameters(14).Fixed = FixorNot(14,1); % r   [m]
nlgrm.parameters(15).Fixed = FixorNot(15,1); % b_r [kg*m^2*s^-1]
nlgrm.parameters(16).Fixed = FixorNot(16,1); % b_l [kg*m^2*s^-1]
nlgrm.parameters(17).Fixed = FixorNot(17,1); % b_p [kg*m^2*s^-1]

%% Visualise model configuration

prompt = 'Do you want visualise extra configurations of the model? Y/N [N]: ';
strvis = input(prompt,'s');
if isempty(strvis)
    strvis = 'N';
end

if strcmp(strvis, 'Y') || strcmp(strvis, 'y')

    % View the inital model
    % a. Get basic information about the model
    % The Otbot has 12 states and 15 model parameters
    disp('This is the dimention of our Otbot model');
    size(nlgrm)

    % b. View the initial states and parameters
    % First we view the initial States
    prompt = 'Do you want to display the inital states configuration? Y/N [N]: ';
    strinput = input(prompt,'s');
    if isempty(strinput)
        strinput = 'N';
    end

    if strcmp(strinput, 'Y') || strcmp(strinput, 'y')
        disp('Intitial state of x')
        nlgrm.InitialStates(1)

        disp('Intitial state of y')
        nlgrm.InitialStates(2)

        disp('Intitial state of alpha')
        nlgrm.InitialStates(3)

        disp('Intitial state of varphi_r')
        nlgrm.InitialStates(4)

        disp('Intitial state of varphi_l')
        nlgrm.InitialStates(5)

        disp('Intitial state of varphi_p')
        nlgrm.InitialStates(6)

        disp('Intitial state of x_dot')
        nlgrm.InitialStates(7)

        disp('Intitial state of y_dot')
        nlgrm.InitialStates(8)

        disp('Intitial state of alpha_dot')
        nlgrm.InitialStates(9)

        disp('Intitial state of varphi_dot_r')
        nlgrm.InitialStates(10)

        disp('Intitial state of varphi_dot_l')
        nlgrm.InitialStates(11)

        disp('Intitial state of varphi_dot_p')
        nlgrm.InitialStates(12)
    elseif strcmp(strinput, 'N') || strcmp(strinput, 'n')
        disp('Omitting the display of the inital states')
    else
        disp('Input not recognised, omitting display of the initial states')
    end

    % Then the parameters 
    prompt = 'Do you want to display the parameters configuration? Y/N [N]: ';
    strinput = input(prompt,'s');
    if isempty(strinput)
        strinput = 'N';
    end

    if strcmp(strinput, 'Y') || strcmp(strinput, 'y')
        nlgrm.Parameters(1)
        nlgrm.Parameters(2)
        nlgrm.Parameters(3)
        nlgrm.Parameters(4)
        nlgrm.Parameters(5)
        nlgrm.Parameters(6)
        nlgrm.Parameters(7)
        nlgrm.Parameters(8)
        nlgrm.Parameters(9)
        nlgrm.Parameters(10)
        nlgrm.Parameters(11)
        nlgrm.Parameters(12)
        nlgrm.Parameters(13)
        nlgrm.Parameters(14)
        nlgrm.Parameters(15)
        nlgrm.Parameters(16)
        nlgrm.Parameters(17)
    elseif strcmp(strinput, 'N') || strcmp(strinput, 'n')
        disp('Omitting display of parameters')
    else
        disp('Input not recognised, omitting display of parameters')
    end

    % Then the initial states that are fixed
    prompt = 'Do you want to display which initial states are Fixed? Y/N [N]: ';
    strinput = input(prompt,'s');
    if isempty(strinput)
        strinput = 'N';
    end

    if strcmp(strinput, 'Y') || strcmp(strinput, 'y')
        getinit(nlgrm,'Fixed')
    elseif strcmp(strinput, 'N') || strcmp(strinput, 'n')
        disp('Omitting display of fixed initial states')
    else
        disp('Input not recognised, omitting display of fixed initial states')
    end


    % Then the minimum values for the parameters
    prompt = 'Do you want to display the minimun allowed values for the parameters? Y/N [N]: ';
    strinput = input(prompt,'s');
    if isempty(strinput)
        strinput = 'N';
    end

    if strcmp(strinput, 'Y') || strcmp(strinput, 'y')
        getpar(nlgrm,'Min')
    elseif strcmp(strinput, 'N') || strcmp(strinput, 'n')
        disp('Omitting display of minimun allowed values for the parameters')
    else
        disp('Input not recognised, omitting display of minimun allowed values for the parameters')
    end

    % Obtain basic information about the object
    disp('Displaying basic information about the idngrey object')
    disp(nlgrm)

    prompt = 'Do you want extra information about the current object? Y/N [N]: ';
    strinput = input(prompt,'s');
    if isempty(strinput)
        strinput = 'N';
    end

    if strcmp(strinput, 'Y') || strcmp(strinput, 'y')
        disp('Running the commands get and present to gain extra information')
        get(nlgrm)
        present(nlgrm)
    elseif strcmp(strinput, 'N') || strcmp(strinput, 'n')
        disp('Omitting display of extra information of the idngrey object')
    else
        disp('Input not recognised, omitting display of extra information of the idngrey object')
    end
elseif strcmp(strvis, 'N') || strcmp(strvis, 'n')
    disp('Omitting the display of any configuration')
else
    disp('Input not recognised, omitting display of any configuration')
end

%% Set up the nlgreyest options
% In this case we only want to choose to display the estimation progress. 
% But there are a lot more options that can be changed as desired. See nlgreyestOptions.
opt = nlgreyestOptions('Display','on');

switch NLGreySearchConfig
    case 0
        disp('Runing nlgreyest with default searching configuration (lsqnonlin)')
    case 1
        disp('Runing nlgreyest with fmincon & interior-point')
        opt.SearchMethod = 'fmincon';
        opt.SearchOptions.Algorithm = 'interior-point';
    otherwise
        disp('WARNING: this NLGreySearchConfig does not exist, runing with the default configuration (lsqnonlin)')
end
% opt.SearchOptions.FunctionTolerance = 1e-30;

%% Estimate the model
% The nlgreyest command updates the parameter of nlgrm.
prompt = 'Do you want to continue with the execution of the script? Y/N [Y]: ';
strinput = input(prompt,'s');
if isempty(strinput)
    strinput = 'Y';
end

if strcmp(strinput, 'Y') || strcmp(strinput, 'y')
    [pvecinit, pvec_sdinit] = getpvec(nlgrm); 
    tStart = tic;
    nlgrm = nlgreyest(data,nlgrm,opt);
    tEnd = toc(tStart); 
    CESFlag = 1;
elseif strcmp(strinput, 'N') || strcmp(strinput, 'n')
    disp('The execution will be terminated')
    CESFlag = 0;
else
    disp('Unrecognised input, terminating the execution')
    CESFlag = 0;
end




