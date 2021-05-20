%% This script is made to save the required data to make the sensitivity plots

filename = 'sensitivity_data';

% Compute the % of deviation
Dev_percentage = com_ini_dev(ParameterRealValues, ParameterInitialGuesses, TSinputDev); 

% Setting up the Free or fixed parameters for estimation vector

FixedVector = [nlgrm.parameters(1).Fixed;   % I [kg*m^2] 
                nlgrm.parameters(2).Fixed]; % b [kg*m^2*s^-1]
            
%% Saveing all the Data in a .mat file

percent = TSinputDev;
percent = num2str(percent,'%.0f');
Tsvalue = num2str(data.Ts,'%.3f');

fullsavename = strcat('SensitivityData/',filename,Tsvalue,'s_',percent,'.mat');

save(fullsavename,"Dev_percentage","FixedVector","RelativeError_Percentage")

%% Extra functions in this document

function Dev_C = com_ini_dev(RVal, IGuess, TSinputDev)

lsv = max(size(RVal));

Dev_C = zeros(lsv,1);

for i=1:lsv
    if RVal(i,1)~=0
        Dev_C(i,1) = abs(RVal(i,1) - IGuess(i,1))./abs(RVal(i,1))*100;
    else
        Dev_C(i,1) = TSinputDev; % We give the value inputed directly
    end
end

end

