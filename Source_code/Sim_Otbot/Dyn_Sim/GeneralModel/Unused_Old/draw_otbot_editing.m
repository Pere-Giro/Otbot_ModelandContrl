function [  ] = draw_otbot_editing( m, q )
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
trplot(eye(3), 'length', 1, 'color', 'k');
hold on;
grid on;

%%%%%%%%%%%%%%%%%
view([0,90]); % this part sets the view of the plot coment if you want isometric

%%%%%%%%%%%%%%%%%
axis([-3 3 -3 3]);

xticks(-100:1:100);
yticks(-100:1:100);

% Drawing chassis Body

% Matrix to define the body in relative frame with the pivot point at the
% origin
CBprel = [0, -m.l_1, -m.l_1;
          0, -m.l_2, +m.l_2;
          zeros(1,3)];
   
% Now we rotate this vectors angle alpha-varphi_p
Rmat = rotz(rad2deg(xs(3)-xs(6))); % Waring robotics toolbox need the angles in degrees
CBp = zeros(3,3);

for i=1:max(length(CBprel))
    CBp(:,i) = Rmat*CBprel(:,i); % Rotation
    CBp(:,i) = CBp(:,i) + [xs(1);xs(2);0]; % Translation
end

fill3(CBp(1,:),CBp(2,:),CBp(3,:),(1/255)*[191,191,191]);

% Drawing Wheels

% Matrix to define the body in relative frame with the pivot point at the
% origin
       
RWBprel = [ + m.r,        + m.r,        - m.r,        - m.r ;
           + (m.l_2)/6, - (m.l_2)/6,  - (m.l_2)/6, + (m.l_2)/6;
           zeros(1,4)];
       
% Now we apply the rotation operation
RWBp_r = zeros(3,4);
RWBp_l = zeros(3,4);
for i=1:max(length(RWBprel))
    RWBp_r(:,i) = Rmat*RWBprel(:,i);  % Rotation
    RWBp_l(:,i) = RWBp_r(:,i);        % Rotation is the same for both wheels 
    RWBp_r(:,i) = RWBp_r(:,i) + CBp(:,2); % Translation for the right wheel
    RWBp_l(:,i) = RWBp_l(:,i) + CBp(:,3); % Translation for the left wheel
end

fill3(RWBp_r(1,:),RWBp_r(2,:),RWBp_r(3,:),(1/255)*[68,68,68]); % Drawing the right wheel
fill3(RWBp_l(1,:),RWBp_l(2,:),RWBp_l(3,:),(1/255)*[68,68,68]); % Drawing the left wheel

% Plot the platform

otbot_circle(xs(1),xs(2),0.3);

% Drawing the orientation line
plot3([xs(1),0.3*cos(xs(3))+xs(1)],[xs(2),0.3*sin(xs(3))+xs(2)],[0,0],'color','r')


hold off;
end

