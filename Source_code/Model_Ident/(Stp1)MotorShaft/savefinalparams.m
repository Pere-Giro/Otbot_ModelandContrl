%% This script is made to save the final values of the estimated parameters

filename = 'final_pars';

percent = TSinputDev;
percent = num2str(percent,'%.0f');
Tsvalue = num2str(data.Ts,'%.3f');

fullsavename = strcat(PartialDirNameExcel,'/',filename,Tsvalue,'s_',percent,'.mat');
fullsavename2 = strcat('Final_Pars_Sim/Motor_Shaft/SetParameters/',filename,'.mat');

save(fullsavename,"ParameterFinalValues")
save(fullsavename2,"ParameterFinalValues")

%% Now we will be saveing the Tsvalue, inputDeviation and Tsaction

Tsvalnum = data.Ts;
save('Final_Pars_Sim/Motor_Shaft/a0_Tsvalnum.mat',"Tsvalnum")
save('Final_Pars_Sim/Motor_Shaft/a1_inputDev.mat',"TSinputDev")

% Save the time of action holding haction or Tsaction
prompt = 'Introduce the time value of holding action, haction or Tsaction [0.01]: ';
Tsaction = input(prompt);
if isempty(Tsaction)
    Tsaction = 0.01;
end

save('Final_Pars_Sim/Motor_Shaft/a3_Tsaction.mat',"Tsaction")

% Save the time of simulation
prompt = 'Introduce total time of simulation [1]: ';
Tsim = input(prompt);
if isempty(Tsim)
    Tsim = 1;
end

save('Final_Pars_Sim/Motor_Shaft/a4_Tsim.mat',"Tsim")

% Saveing Sampling instants
SampInst = data.SamplingInstants;
save('Final_Pars_Sim/Motor_Shaft/a5_SampInstants.mat',"SampInst")
