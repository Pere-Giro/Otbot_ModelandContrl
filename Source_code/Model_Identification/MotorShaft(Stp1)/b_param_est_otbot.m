%% Specify Initial Sate estimation
% By default the inital states will not be estimated if we want to estimate
% them we must set it up.

nlgrm.InitialStates(1).Fixed = true;  % varphi          [rad]
nlgrm.InitialStates(2).Fixed = true;  % varphi_dot      [rad/s]

%% Specify known parameters
% Specify which parameters are alredy known and which wants we need to
% estimate

nlgrm.parameters(1).Fixed = FixorNot(1,1);  % I [kg*m^2]
nlgrm.parameters(2).Fixed = FixorNot(2,1);  % b [kg*m^2*s^-1]

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
% opt.SearchOptions.FunctionTolerance = 1e-10;

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
    nlgrm = nlgreyest(data,nlgrm, opt);
    tEnd = toc(tStart); 
    CESFlag = 1;
elseif strcmp(strinput, 'N') || strcmp(strinput, 'n')
    disp('The execution will be terminated')
    CESFlag = 0;
else
    disp('Unrecognised input, terminating the execution')
    CESFlag = 0;
end




