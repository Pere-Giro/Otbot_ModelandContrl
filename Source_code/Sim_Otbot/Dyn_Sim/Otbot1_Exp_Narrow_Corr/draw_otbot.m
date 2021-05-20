function [  ] = draw_otbot( m, q )
%DRAW_DM Summary of this function goes here
%   Detailed explanation goes here
figure(100);
pos_fig1 = [0 0 1920 1080];
set(gcf,'Position',pos_fig1)

xs(1)=q.x;
xs(2)=q.y;
xs(3)=q.alpha;
xs(4)=q.varphi_r;
xs(5)=q.varphi_l;
xs(6)=q.varphi_p;

% Drawing axis
% trplot(eye(3), 'length', 1, 'color', 'k');
clf(figure(100))
hold on;
grid on;
box on;
axis equal

%%%%%%%%%%%%%%%%%
view([0,90]); % this part sets the view of the plot coment if you want isometric

%%%%%%%%%%%%%%%%%
axis([-1 10 -1 4]);

xticks(-100:1:100);
yticks(-100:1:100);

% Defining rotation matrix
Rmat = m.Rzmat(xs(3)-xs(6));

Rmat2 = m.Rzmat(xs(3));

% Drawing chassis Body

trvec = [xs(1);xs(2);0];

CBp = Rmat*m.CBprel + [trvec,trvec,trvec]; % Rotation + Translation

fill3(CBp(1,:),CBp(2,:),CBp(3,:),(1/255)*[191,191,191]);

% Drawing Wheels     
RWBp_r = Rmat*m.RWBprel; % Rotation

RWBp_l = RWBp_r;   % Rotation is the same for both wheels  

RWBp_r = RWBp_r + [CBp(:,2),CBp(:,2),CBp(:,2),CBp(:,2)]; % Translation 
RWBp_l = RWBp_l + [CBp(:,3),CBp(:,3),CBp(:,3),CBp(:,3)]; % Translation

fill3(RWBp_r(1,:),RWBp_r(2,:),RWBp_r(3,:),(1/255)*[68,68,68]); % Drawing the right wheel
fill3(RWBp_l(1,:),RWBp_l(2,:),RWBp_l(3,:),(1/255)*[68,68,68]); % Drawing the left wheel

% Plot the platform

otbot_circle_v2(xs(1),xs(2),0.45,'NO');

%----------Draw the path-------------%

% Piece1
plot([0,3],[0,0],'g', 'LineWidth', 1)

% Piece2
plot([3,3],[0,3],'g', 'LineWidth', 1)

% Piece3
plot([3,6],[3,3],'g', 'LineWidth', 1)

% Piece4
plot([6,6],[3,0],'g', 'LineWidth', 1)

% Piece5
plot([6,9],[0,0],'g', 'LineWidth', 1)

% Drawing the orientation line
plot3([xs(1),0.45*cos(xs(3))+xs(1)],[xs(2),0.45*sin(xs(3))+xs(2)],[0,0],'color','r')

trih = 0.45/6;
triedge = trih/cos(pi/6);
tricords0 = [0,             trih  ,               trih;
             0, triedge*sin(pi/6) , -triedge*sin(pi/6);
             zeros(1,3)]; % Triangle coords at the origin
         
trirot = Rmat2*tricords0; % Rotate triangle
tricords = trirot + [(0.45-trih)*cos(xs(3))+xs(1);(0.45-trih)*sin(xs(3))+xs(2);0];
fill(tricords(1,:),tricords(2,:),'r','EdgeColor','none') % Drawing

%---------------------------------------------------%
% Plotting center of mass of the chassis body
% GBp = Rmat*[m.x_G; m.y_G; 0] + trvec;
% 
% otbot_circle_v2(GBp(1),GBp(2),m.l_2/6,'b');
%---------------------------------------------------%

% Plotting center of mass of the platform
GPp = Rmat2*[m.x_F; m.y_F; 0] + trvec;

otbot_circle_v2(GPp(1),GPp(2),m.l_2/6,'g');

%----------Draw the walls-------------%

% Define the color
% pale_grey = [125 125 125]/255; % defined color but not used in the end we will
                               % be using red

% Define the walls
poly1 = [-2,   -2,   11, 11;
         -2, -0.6, -0.6, -2];
     
poly2 = [ 3.6, 3.6, 5.4,  5.4;
         -0.6, 2.4, 2.4, -0.6];

poly3 = [ -2, -2, 2.4, 2.4;
         0.6,  5,   5, 0.6];

poly4 = [6.6, 6.6, 11,  11;
         0.6,   5,  5, 0.6];

poly5 = [ -2, -2, 11,  11;
         3.6,  5,  5, 3.6];
     
% Drawing
     
fill(poly1(1,:),poly1(2,:),'r','EdgeColor','None')

fill(poly2(1,:),poly2(2,:),'r','EdgeColor','None')

fill(poly3(1,:),poly3(2,:),'r','EdgeColor','None')

fill(poly4(1,:),poly4(2,:),'r','EdgeColor','None')

fill(poly5(1,:),poly5(2,:),'r','EdgeColor','None')

rectangle('Position', [-1,-1, 11,5])

%----------Draw the load (Filled Square) -------------%

% % Offets form the center of pivot
% ofsetval = 0.42;
% 
% % Color
% pale_yellow = [255 250 151]/255;
% 
% % Coordinates of the load
% 
% LoaCor = [GPp(1)-ofsetval, GPp(1)-ofsetval, GPp(1)+ofsetval, GPp(1)+ofsetval;
%           GPp(2)-ofsetval, GPp(2)+ofsetval, GPp(2)+ofsetval, GPp(2)-ofsetval];
%            
% fill(LoaCor(1,:),LoaCor(2,:),pale_yellow)

%----------Draw the load (X-ray)-------------%

% Offets form the center of pivot
ofsetval = 0.42;

% Coords for the rectangle x,y,w,h 

xywh = [GPp(1)-ofsetval,GPp(2)-ofsetval,2*ofsetval,2*ofsetval];

rectangle('Position',xywh)
     
hold off;
end

