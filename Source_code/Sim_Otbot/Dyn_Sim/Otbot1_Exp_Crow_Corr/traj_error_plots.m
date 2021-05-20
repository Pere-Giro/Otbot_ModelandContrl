close all
load('Crow_Corr_Data/OLstatesCOMpp')

%% Set the Flags and the name

CN=' Crowded_Corridor_Test2_01OL'; % Name the image (Warning: avoid dots in the name if EPS save)

Flag_SEP = 'YES'; % Flag to choose if one wants to save the plots or not
                 % (Flag_SEP == YES)> We save the plots
                 % (Flag_SEP == NO)> We do not save the plots
                 
Flag_EPS = 'YES'; % Flag to choose if the images will be saves in eps or png)
                % (Flag_EPS == 'YES' > Images saved in EPS)
                % (Flag_EPS == 'NO'  > Images saved in png)


%% Compute the errors of trajectory in each of the task coordinates

x_error = OLstatesCoMpp(1,:) - sstates(1,1:9991);
y_error = OLstatesCoMpp(2,:) - sstates(2,1:9991);
alpha_error = OLstatesCoMpp(3,:) - sstates(3,1:9991);

xdot_error = OLstatesCoMpp(7,:) - sstates(4,1:9991);
ydot_error = OLstatesCoMpp(8,:) - sstates(5,1:9991);
alphadot_error = OLstatesCoMpp(9,:) - sstates(6,1:9991);

% Plot the errors verison Full rectangle

% figure('WindowState','maximized');
% plot(times,x_error)
% xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$x$ error [m]','Interpreter','latex','FontSize', 17),title('Configuration error','Interpreter','latex','FontSize', 17)
% grid on
% 
% figure('WindowState','maximized');
% plot(times,y_error)
% xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$y$ error [m]','Interpreter','latex','FontSize', 17),title('Configuration error','Interpreter','latex','FontSize', 17)
% grid on
% 
% figure('WindowState','maximized');
% plot(times,alpha_error)
% xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\alpha$ error [rad]','Interpreter','latex','FontSize', 17),title('Configuration error','Interpreter','latex','FontSize', 17)
% grid on
% 
% figure('WindowState','maximized');
% plot(times,xdot_error)
% xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\dot{x}$ error [m/s]','Interpreter','latex','FontSize', 17),title('Velocity error','Interpreter','latex','FontSize', 17)
% grid on
% 
% figure('WindowState','maximized');
% plot(times,ydot_error)
% xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\dot{y}$ error [m/s]','Interpreter','latex','FontSize', 17),title('Velocity error','Interpreter','latex','FontSize', 17)
% grid on
% 
% figure('WindowState','maximized');
% plot(times,alphadot_error)
% xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('$\dot{\alpha}$ error [rad/s]','Interpreter','latex','FontSize', 17),title('Velocity error','Interpreter','latex','FontSize', 17)
% grid on
% 
% % Plot error x and y in same plot
% figure('WindowState','maximized');
% plot(times,x_error)
% hold on
% plot(times,y_error)
% hold off
% xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('error [m]','Interpreter','latex','FontSize', 17),title('Configuration Error','Interpreter','latex','FontSize', 17)
% legend('$x$ error','$y$ error','Interpreter','latex','FontSize', 17)
% grid on
% 
% % Plot error x_dot and y_dot in same plot
% figure('WindowState','maximized');
% plot(times,xdot_error)
% hold on
% plot(times,ydot_error)
% hold off
% xlabel('t [s]','Interpreter','latex','FontSize', 17),ylabel('error [m/s]','Interpreter','latex','FontSize', 17),title('Velocity Error','Interpreter','latex','FontSize', 17)
% legend('$\dot{x}$ error','$\dot{y}$ error','Interpreter','latex','FontSize', 17)
% grid on


% Plot the errors verison small square

figure;
plot(times,x_error)
xlabel('t [s]','Interpreter','latex' ),ylabel('$x$ error [m]','Interpreter','latex' ),title('Configuration error','Interpreter','latex' )
grid on

figure;
plot(times,y_error)
xlabel('t [s]','Interpreter','latex' ),ylabel('$y$ error [m]','Interpreter','latex' ),title('Configuration error','Interpreter','latex' )
grid on

figure;
plot(times,alpha_error)
xlabel('t [s]','Interpreter','latex' ),ylabel('$\alpha$ error [rad]','Interpreter','latex' ),title('Configuration error in $\alpha$','Interpreter','latex' )
grid on

figure;
plot(times,xdot_error)
xlabel('t [s]','Interpreter','latex' ),ylabel('$\dot{x}$ error [m/s]','Interpreter','latex' ),title('Velocity error','Interpreter','latex' )
grid on

figure;
plot(times,ydot_error)
xlabel('t [s]','Interpreter','latex' ),ylabel('$\dot{y}$ error [m/s]','Interpreter','latex' ),title('Velocity error','Interpreter','latex' )
grid on

figure;
plot(times,alphadot_error)
xlabel('t [s]','Interpreter','latex' ),ylabel('$\dot{\alpha}$ error [rad/s]','Interpreter','latex' ),title('Velocity error in $\dot{\alpha}$','Interpreter','latex' )
grid on

% Plot error x and y in same plot
figure;
plot(times,x_error)
hold on
plot(times,y_error)
hold off
xlabel('t [s]','Interpreter','latex' ),ylabel('error [m]','Interpreter','latex' ),title('Configuration error in $x$ and $y$','Interpreter','latex' )
legend('$x$ error','$y$ error','Interpreter','latex' )
grid on

% Plot error x_dot and y_dot in same plot
figure;
plot(times,xdot_error)
hold on
plot(times,ydot_error)
hold off
xlabel('t [s]','Interpreter','latex' ),ylabel('error [m/s]','Interpreter','latex' ),title('Velocity error in $\dot{x}$ and $\dot{y}$','Interpreter','latex' )
legend('$\dot{x}$ error','$\dot{y}$ error','Interpreter','latex' )
grid on

%% Save the plots

switch Flag_SEP
    case 'YES'
        for i=1:8
            if i==1
                names='15 - x_error';
            elseif i==2
                names='16 - y_error';
            elseif i==3
                names='17 - alpha_error';
            elseif i==4
                names='18 - xdot_error';
            elseif i==5
                names='19 - ydot_error';
            elseif i==6
                names='20 - alphadot_error';
            elseif i==7
                names='21 - Error_xy';
            elseif i==8
                names='22 - Error_dxdy';
            end
            switch Flag_EPS
                case 'YES'
                    saveas(figure(i),strcat(names,CN,'.eps'),'epsc');
                case 'NO'
                    saveas(figure(i),strcat(names,CN,'.png'));
                otherwise
                    disp('Error this Flag_EPS does not exist')
            end
        end
    case 'NO'
        disp('Trajectory error plots will not be saved')
    otherwise
        disp('WARNING: This Flag_SEP does not exist not saveing plots')
end

