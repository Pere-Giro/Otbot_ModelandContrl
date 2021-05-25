%% Step 1 Load data and set-up the idngrey
% we call the script to do so
a_setup_otbot_gbme
%% Step 2 - Set up and Model Identification
% In this section we set up all the free parameters the fixed ones and same
% for inicial states. Also we display all the configurations for the user
% in order to check if everything has been set correctly.
b_param_est_otbot

%% Displaying the results
if CESFlag == 1
    [pvec, pvec_sd] = getpvec(nlgrm);
    load('real_pars_data.mat')
    ParameterRealValues = real_pars;
    ParameterInitialGuesses = pvecinit;
    ParameterFinalValues = pvec;
    
    prompt = 'Estimation terminated successfully! Do you want to display the results? Y/N [Y]: ';
    mainstrinput = input(prompt,'s');
    if isempty(mainstrinput)
        mainstrinput = 'Y';
    end
    
    ParameterName = {'I_b'; 'I_p'; 'I_a'; 'I_t'; 'l_1'; 'l_2'; 'm_b'; 'm_w'; 'm_p'; 'x_G'; 'y_G'; 'x_F'; 'y_F'; 'r'; 'Abias'; 'Xbias'; 'Ybias'};
%         RelativeError_Percentage = abs(ParameterRealValues - ParameterFinalValues)./ParameterRealValues*100;
    [RelativeError_Percentage, AbsoluteError] = rel_error(ParameterRealValues, ParameterFinalValues);
    StDev = pvec_sd;
    Units = {'kg*m^2'; 'kg*m^2'; 'kg*m^2'; 'kg*m^2'; 'm'; 'm'; 'kg'; 'kg'; 'kg'; 'm'; 'm'; 'm'; 'm'; 'm'; 'nounit'; 'nounit'; 'nounit'};
    ResultsTable = table(ParameterName, Units, ParameterRealValues, ParameterInitialGuesses, ParameterFinalValues, AbsoluteError, RelativeError_Percentage, StDev );
        
    if strcmp(mainstrinput, 'Y') || strcmp(mainstrinput, 'y')
        disp('Displaying a Table with the final values of the parameters')
        disp(ResultsTable)
    elseif strcmp(mainstrinput, 'N') || strcmp(mainstrinput, 'n')
        disp('The display of results will be omitted')
    else
        disp('Unrecognised input, omitting the display of results')
    end
    
    % Save everything?
    prompt = 'Do you want to save everything? Y/N [Y]: ';
    mainstrinput = input(prompt,'s');
    if isempty(mainstrinput)
        mainstrinput = 'Y';
    end

    if strcmp(mainstrinput, 'Y') || strcmp(mainstrinput, 'y')
        savefinalparams
        save_exceltable
%         writetable(ResultsTable,'grey_est_results.xlsx')
        save_sensitivitydata
        save_Report
    elseif strcmp(mainstrinput, 'N') || strcmp(mainstrinput, 'n')
        % Save the final values of the parameters
        prompt = 'Do you want to save final values of the Parameters in .mat file? Y/N [N]: ';
        mainstrinput = input(prompt,'s');
        if isempty(mainstrinput)
            mainstrinput = 'N';
        end

        if strcmp(mainstrinput, 'Y') || strcmp(mainstrinput, 'y')
            savefinalparams
        elseif strcmp(mainstrinput, 'N') || strcmp(mainstrinput, 'n')
            disp('The table of results will not be saved')
        else
            disp('Unrecognised input, omitting the display of results')
        end

        % Save the table in Excel

        prompt = 'Do you want to save the results table in Excel? Y/N [N]: ';
        mainstrinput = input(prompt,'s');
        if isempty(mainstrinput)
            mainstrinput = 'N';
        end

        if strcmp(mainstrinput, 'Y') || strcmp(mainstrinput, 'y')
            save_exceltable
%         writetable(ResultsTable,'grey_est_results.xlsx')
        elseif strcmp(mainstrinput, 'N') || strcmp(mainstrinput, 'n')
            disp('The table of results will not be saved')
        else
            disp('Unrecognised input, omitting the display of results')
        end

        % Save the data for sensitivity analisis

        prompt = 'Do you want to save the data for the sensitivity analisis? Y/N [N]: ';
        mainstrinput = input(prompt,'s');
        if isempty(mainstrinput)
            mainstrinput = 'N';
        end

        if strcmp(mainstrinput, 'Y') || strcmp(mainstrinput, 'y')
            save_sensitivitydata
        elseif strcmp(mainstrinput, 'N') || strcmp(mainstrinput, 'n')
            disp('The data for sensitivity analisis will not be saved')
        else
            disp('Unrecognised input, data for sensitivity analisis will not be saved')
        end

        % Save Results Report

        prompt = 'Do you want to save results report? Y/N [N]: ';
        mainstrinput = input(prompt,'s');
        if isempty(mainstrinput)
            mainstrinput = 'N';
        end

        if strcmp(mainstrinput, 'Y') || strcmp(mainstrinput, 'y')
            save_Report
        elseif strcmp(mainstrinput, 'N') || strcmp(mainstrinput, 'n')
            disp('The results report will not be saved')
        else
            disp('Unrecognised input, results report will not be saved')
        end
    else
       disp('Unrecognised input, nothing will be saved')
    end
end