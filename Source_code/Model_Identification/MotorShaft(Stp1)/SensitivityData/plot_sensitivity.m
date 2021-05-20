%% Loading all the .mat files

SPN='SP_Data0.010s_Test1.01_05s.png'; % Name the image

Tsplot = '0.010';
FilenamePlot = strcat('Sensitivity_Data',Tsplot,'s_'); 
DataPlot(2,3,1) = 0;
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
ParNames = ["Inertia of the motor shaft [kg$\cdot$m$^2$]";
            "Viscous friction coefficient [kg$\cdot$m$^3 \cdot s^{-1}$]"];

%% Rearrange the arrays create auxiliar variables
Dim = size(DataPlot);
% Compute the numer of free parameters
NotFix = ~DataPlot(:,2,1);
NumPar = ones(Dim(1),1);
NumPar = NumPar(NotFix);
NumPar = sum(NumPar);
ParNamesList = ParNames(NotFix);

% Creating the final parameter array to plot
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



