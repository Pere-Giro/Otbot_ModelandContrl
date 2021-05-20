%% Loading all the .mat files

SPN='SP_Data0.010s_Test1.02_3s.png'; % Name the image

Tsplot = '0.010';
FilenamePlot = strcat('Sensitivity_Data',Tsplot,'s_'); 
DataPlot(17,3,3) = 0;
jnum = 1;
for ipercent=10:10:50
    Percentvaluename = num2str(ipercent,'%.0f');
    FullFilename = strcat(FilenamePlot,Percentvaluename,'.mat');
    load(FullFilename);
    DataPlot(:,1,jnum) = Dev_percentage;
    DataPlot(:,2,jnum) = FixedVector;
    DataPlot(:,3,jnum) = RelativeError_Percentage;
    jnum = jnum+1;
end

% Create parameter names list
ParNames = ["Inertia of the chassis [kg$\cdot$m$^2$]";
            "Inertia of the platform [kg$\cdot$m$^2$]";
            "Axial inertia of one wheel [kg$\cdot$m$^2$]";
            "Twisting inertia of one wheel [kg$\cdot$m$^2$]";
            "Pivot offset [m]";
            "One half of the wheels separation [m]";
            "Mass of the chassis [kg]";
            "Mass of one wheel [kg]";
            "Mass of the platform [kg]";
            "x coord of chassis [m]";
            "y coord of chassis [m]";
            "x coord of platform body [m]";
            "y coord of platform body [m]";
            "Wheel radius [m]";
            "Friction coef. right wheel [kg$\cdot$m$^2\cdot s^{-1}$]";
            "Friction coef. left wheel [kg$\cdot$m$^2\cdot s^{-1}$]";
            "Friction coef. pivot joint [kg$\cdot$m$^2\cdot s^{-1}$]"];

%% Rearrange the arrays create auxiliar variables
Dim = size(DataPlot);
% Compute the numer of free parameters
NotFix = ~DataPlot(:,2,1);
NumPar = ones(Dim(1),1);
NumPar = NumPar(NotFix);
NumPar = sum(NumPar);
ParNamesList = ParNames(NotFix);

%Creating the final parameter array to plot
Params(2,Dim(3),NumPar) = 0;
count1 = 1;
for i=1:Dim(3)
    count = 1;
    for j=1:Dim(1)
        if ~DataPlot(j,2,i)
            Params(1,count1,count) = DataPlot(j,1,i);
            Params(2,count1,count) = DataPlot(j,3,i);
            count = count + 1;   
        end
    end 
    count1 = count1 + 1;
end

%% Plotting
f = figure('WindowState','maximized');
legend('Interpreter','latex')
grid on
hold on
xticks(-100:10:100);
for i=1:NumPar
    plot(Params(1,:,i),Params(2,:,i),'o-','DisplayName',ParNamesList(i))
end
xlabel('Deviation of initial guess [\%]','Interpreter','latex'),ylabel('Relative Error [\%]','Interpreter','latex'),title('Sensitivity Graphic','Interpreter','latex')
hold off

%% Saveing the plot
saveas(gcf,SPN);



