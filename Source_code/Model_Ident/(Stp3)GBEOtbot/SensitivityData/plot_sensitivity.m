%% Loading all the .mat files

SPN=' SP_Data0_010s_Test1_01_1s'; % Name the image

Flag_EPS = 'YES'; % Flag to choose if the images will be saves in eps or png)
                % (Flag_EPS == 'YES' > Images saved in EPS)
                % (Flag_EPS == 'NO'  > Images saved in png)

Tsplot = '0.010';
FilenamePlot = strcat('Sensitivity_Data',Tsplot,'s_'); 
% With this we manage the .mat files that will be loaded
sizevec = 1:1:32;

sizeDP = size(sizevec,2);
DataPlot(17,3,sizeDP) = 0;
jnum = 1;
for ipercent = sizevec
    Percentvaluename = num2str(ipercent,'%.0f');
    FullFilename = strcat(FilenamePlot,Percentvaluename,'.mat');
    load(FullFilename);
    DataPlot(:,1,jnum) = Dev_value;
    DataPlot(:,2,jnum) = FixedVector;
    DataPlot(:,3,jnum) = RelativeError_Percentage;
    DataPlot(:,4,jnum) = AbsoluteError;
    jnum = jnum+1;
end

% Create parameter names list [This list is for the first 2 ways of doing
% the plots with legend]

% ParNames = ["Inertia of the chassis body [kg$\cdot$m$^2$]";
%             "Inertia of the platform [kg$\cdot$m$^2$]";
%             "Axial inertia of one wheel [kg$\cdot$m$^2$]";
%             "Twisting inertia of one wheel [kg$\cdot$m$^2$]";
%             "Pivot offset [m]";
%             "One half of the wheels separation [m]";
%             "Mass of the chassis base [kg]";
%             "Mass of one wheel [kg]";
%             "Mass of the platform [kg]";
%             "x coord of chassis body [m]";
%             "y coord of chassis body [m]";
%             "x coord of platform body [m]";
%             "y coord of platform body [m]";
%             "Wheel radius [m]";
%             "Alphadot Bias [nounit]";
%             "xddot Bias [nounit]";
%             "yddot Bias [nounit]"];

% Create parameter names list [This list is for the 3rd unlimate way of doing
% the plots with no legend]

ParNames = ["Inertia of the chassis body";
            "Inertia of the platform";
            "Axial inertia of one wheel";
            "Twisting inertia of one wheel";
            "Pivot offset";
            "One half of the wheels separation";
            "Mass of the chassis base";
            "Mass of one wheel";
            "Mass of the platform";
            "x coord. of the C.O.M. of the chassis body";
            "y coord. of the C.O.M. of chassis body";
            "x coord. of the C.O.M. of platform";
            "y coord. of the C.O.M. of platform";
            "Wheel radius";
            "Alphadot Bias";
            "xddot Bias";
            "yddot Bias"];
        
ParNamesShort = ["I_b";
                 "I_p";
                 "I_a";
                 "I_t";
                 "l_1";
                 "l_2";
                 "m_b";
                 "m_w";
                 "m_p";
                 "x_G";
                 "y_G";
                 "x_F";
                 "y_F";
                 "r";
                 "Abias";
                 "Xbias";
                 "Ybias"];
             
UnitsVect = ["[kg$\cdot$m$^2$]";
             "[kg$\cdot$m$^2$]";
             "[kg$\cdot$m$^2$]";
             "[kg$\cdot$m$^2$]";
             "[m]";
             "[m]";
             "[kg]";
             "[kg]";
             "[kg]";
             "[m]";
             "[m]";
             "[m]";
             "[m]";
             "[m]";
             "[nounit]";
             "[nounit]";
             "[nounit]"];

%% Rearrange the arrays create auxiliar variables
Dim = size(DataPlot);
% Compute the numer of free parameters
NotFix = ~DataPlot(:,2,1);
NumPar = ones(Dim(1),1);
NumPar = NumPar(NotFix);
NumPar = sum(NumPar);
ParNamesList = ParNames(NotFix);
UnitNameList = UnitsVect(NotFix);

% Creating the final parameter array to plot
Params(3,Dim(3),NumPar) = 0;
count1 = 1;
for i=1:Dim(3)
    count = 1;
    for j=1:Dim(1)
        if ~DataPlot(j,2,i)
            Params(1,count1,count) = DataPlot(j,1,i);
            Params(2,count1,count) = DataPlot(j,3,i);
            Params(3,count1,count) = DataPlot(j,4,i);
            count = count + 1;   
        end
    end 
    count1 = count1 + 1;
end

%% Plotting

% old way 'classic' to do the plots
% for i=1:NumPar
%     figure('WindowState','maximized');
%     plot(Params(1,:,i),Params(3,:,i),'o-','DisplayName',ParNamesList(i))
%     xlabel('Deviation of initial guess','Interpreter','latex'),ylabel('Absolute error','Interpreter','latex'),title('Sensitivity Graphic','Interpreter','latex')
%     legend('Interpreter','latex')
%     grid on
%     box on
% %     xticks(-100:10:100);
% end

% New way is more cute for the TFM
% for i=1:NumPar
%     figure('WindowState','maximized');
%     plot(Params(1,:,i),Params(3,:,i),'o-','DisplayName',ParNamesList(i))
%     xlabel('Deviation of initial guess from the real value','Interpreter','latex','FontSize', 17),ylabel('Absolute error','Interpreter','latex','FontSize', 17),title('Sensitivity plot','Interpreter','latex','FontSize', 17)
%     legend('Interpreter','latex','FontSize', 17)
%     grid on
%     box on
% end

% Ultimate way of doing the plots for the TFM
for i=1:NumPar
    figure('WindowState','maximized');
    plot(Params(1,:,i),Params(3,:,i),'o-')
    xlabel(strcat("Deviation of initial guess from the real value"," ",UnitNameList(i)),'Interpreter','latex','FontSize', 17),ylabel(strcat("Absolute error"," ",UnitNameList(i)),'Interpreter','latex','FontSize', 17),title(ParNamesList(i),'Interpreter','latex','FontSize', 17)
    % legend('Interpreter','latex','FontSize', 17)
    grid on
    box on
end

%% Saveing the plot
% saveas(gcf,SPN);

PNSL = ParNamesShort(NotFix);
for i=1:NumPar
    idx = num2str(i);
    name = strcat(idx," - ",PNSL(i),SPN);
%     saveas(figure(i),name);
    switch Flag_EPS
        case 'YES'
            saveas(figure(i),strcat(name,'.eps'),'epsc');
        case 'NO'
            saveas(figure(i),strcat(name,'.png'));
        otherwise
            disp('Error this Flag_EPS does not exist')
    end
end


