%% Lets try the arrow function

figure(100);
% pos_fig1 = [0 0 1920 1080];
% pos_fig1 = [0 0 900 900];
% set(gcf,'Position',pos_fig1)

clf(figure(100))
hold on;
grid on;
box on;
axis equal

%%%%%%%%%%%%%%%%%
view([0,90]); % this part sets the view of the plot coment if you want isometric

%%%%%%%%%%%%%%%%%
axis([-2 4 -3 3]);

xticks(-100:1:100);
yticks(-100:1:100);


% Draw the arrow

arrows(0,0,2,0,'FaceColor','r','EdgeColor','none')



