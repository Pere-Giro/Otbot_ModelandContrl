%% Plotting the outputs measured and action inputs of the Data File.
% This are the data that we will be feeding to nlgreyest, so it may be
% interesting to plot in order to see exacly what are how are we building
% up the Constrained Minimization Problem.

%------- Plot the outputs -------%

% Plot alphadot coordinate
figure;
plot(data.SamplingInstants, data.OutputData(:,1),'.')
xlabel('t(s)','Interpreter','latex'),ylabel('$\dot{\alpha}$(rad/s)','Interpreter','latex'),title('Velocity','Interpreter','latex')

% Plot x_ddot coordinate
figure;
plot(data.SamplingInstants, data.OutputData(:,2),'.')
xlabel('t(s)','Interpreter','latex'),ylabel('$\ddot{x}$(m/$s^2$)','Interpreter','latex'),title('Acceleration','Interpreter','latex')

% Plot varphidot_l coordinate
figure;
plot(data.SamplingInstants, data.OutputData(:,3),'.')
xlabel('t(s)','Interpreter','latex'),ylabel('$\ddot{y}$(m/$s^2$)','Interpreter','latex'),title('Acceleration','Interpreter','latex')


%------- Plot the inputs -------%

figure;
plot(data.SamplingInstants,data.InputData(:,1:2),'.')
xlabel('t(s)','Interpreter','latex'),ylabel('$\tau$ of wheels(N$\cdot$m)','Interpreter','latex'),title('Action Torques','Interpreter','latex')
legend('$\tau$ right','$\tau$ left','Interpreter','latex')

figure;
plot(data.SamplingInstants,data.InputData(:,3),'.')
xlabel('t(s)','Interpreter','latex'),ylabel('$\tau$ pivot(N$\cdot$m)','Interpreter','latex'),title('Action Torques','Interpreter','latex')


