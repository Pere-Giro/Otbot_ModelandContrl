%% Star by setting everything up

figure(100);
pos_fig1 = [0 0 1920 1080];
set(gcf,'Position',pos_fig1)

hold on;
grid on;
box on;

axis equal

view([0,90]); % this part sets the view of the plot coment if you want isometric

%%%%%%%%%%%%%%%%%
axis([-1 10 -1 4]);

xticks(-100:1:100);
yticks(-100:1:100);

%% Defining polygons

poly1 = [-2,   -2,   11, 11;
         -2, -0.5, -0.5, -2];
     
poly2 = [ 3.5, 3.5, 5.5,  5.5;
         -0.5, 2.5, 2.5, -0.5];

poly3 = [ -2, -2, 2.5, 2.5;
         0.5,  5,   5, 0.5];

poly4 = [6.5, 6.5, 11,  11;
         0.5,   5,  5, 0.5];

poly5 = [ -2, -2, 11,  11;
         3.5,  5,  5, 3.5];
     
%% Define colors

pale_yellow = [255 250 151]/255;

pale_grey = [125 125 125]/255;

%% Plotting
     
hold on

% Piece1
plot([0,3],[0,0],'g')

% Piece2
plot([3,3],[0,3],'g')

% Piece3
plot([3,6],[3,3],'g')

% Piece4
plot([6,6],[3,0],'g')

% Piece5
plot([6,9],[0,0],'g')

fill(poly1(1,:),poly1(2,:),pale_grey,'EdgeColor','None')

fill(poly2(1,:),poly2(2,:),pale_grey,'EdgeColor','None')

fill(poly3(1,:),poly3(2,:),pale_grey,'EdgeColor','None')

fill(poly4(1,:),poly4(2,:),pale_grey,'EdgeColor','None')

fill(poly5(1,:),poly5(2,:),pale_grey,'EdgeColor','None')

rectangle('Position', [-1,-1, 11,5])


hold off
