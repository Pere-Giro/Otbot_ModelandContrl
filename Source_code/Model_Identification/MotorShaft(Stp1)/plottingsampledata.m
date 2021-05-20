%% Plotting the outputs measured and action inputs of the Data File.
% This are the data that we will be feeding to nlgreyest, so it may be
% interesting to plot in order to see exacly what are how are we building
% up the Constrained Minimization Problem.

%------- Plot the outputs -------%

% Plot varphi_dot coordinate
figure;
plot(data.SamplingInstants, data.OutputData(:,1),'.')
xlabel('t(s)','Interpreter','latex'),ylabel('$\dot{\varphi}$(rad/s)','Interpreter','latex'),title('Velocity','Interpreter','latex')

%------- Plot the input -------%
figure;
plot(data.SamplingInstants,data.InputData(:,1),'.')
xlabel('t(s)','Interpreter','latex'),ylabel('$\tau$ (N$\cdot$m)','Interpreter','latex'),title('Action Torque','Interpreter','latex')


