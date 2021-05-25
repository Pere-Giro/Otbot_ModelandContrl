function [xdot,yout] = otbot1_m_V2_IMU(t, xs, u, I_b, I_p, I_a, I_t, l_1, l_2, m_b, m_w, m_p, x_G, y_G, x_F, y_F, r, b_r, b_l, b_p, varargin)
%OTBOT_M Summary of this function goes here
%  Function containing the equations of motion of Otbot

%------------ List of state variables ------------%
% Note how some of this state variables (the ones commented out) are not 
% used to compute the equations of motion of the system.

% x = xs(1);
% y = xs(2);
alpha = xs(3);
% varphi_r = xs(4);
% varphi_l = xs(5);
varphi_p = xs(6);

% x_dot = xs(7);
% y_dot = xs(8);
alpha_dot = xs(9);
varphi_dot_r = xs(10);
varphi_dot_l = xs(11);
varphi_dot_p = xs(12);

%------------ Create viscous friction term ------------%

tau_friction(1:3,1) = zeros(3,1);
        tau_friction(4,1) = -b_r*varphi_dot_r; % tau_friction_right
        tau_friction(5,1) = -b_l*varphi_dot_l; % tau_friction_left
        tau_friction(6,1) = -b_p*varphi_dot_p; % tau_friction_pivot

%------------ Create system matrices ------------%

MIIKs = reshape([(l_1.*cos(alpha-varphi_p)-l_2.*sin(alpha-varphi_p))./(l_1.*r),(l_1.*cos(alpha-varphi_p)+l_2.*sin(alpha-varphi_p))./(l_1.*r),sin(alpha-varphi_p)./l_1,(l_2.*cos(alpha-varphi_p)+l_1.*sin(alpha-varphi_p))./(l_1.*r),-(l_2.*cos(alpha-varphi_p)-l_1.*sin(alpha-varphi_p))./(l_1.*r),-cos(alpha-varphi_p)./l_1,0.0,0.0,1.0],[3,3]);
M_bars2 = reshape([((I_a.*l_2.^2.*sin(alpha-varphi_p).*2.0+I_b.*r.^2.*sin(alpha-varphi_p)+I_t.*r.^2.*sin(alpha-varphi_p).*2.0+l_1.^2.*m_b.*r.^2.*sin(alpha-varphi_p)+l_1.^2.*m_p.*r.^2.*sin(alpha-varphi_p)+l_2.^2.*m_w.*r.^2.*sin(alpha-varphi_p).*2.0-I_a.*l_1.*l_2.*cos(alpha-varphi_p).*2.0+m_b.*r.^2.*x_G.^2.*sin(alpha-varphi_p)+m_b.*r.^2.*y_G.^2.*sin(alpha-varphi_p)+l_1.*m_p.*r.^2.*y_F.*cos(alpha)+l_1.*m_p.*r.^2.*x_F.*sin(alpha)-l_1.*l_2.*m_b.*r.^2.*cos(alpha-varphi_p)-l_1.*l_2.*m_p.*r.^2.*cos(alpha-varphi_p)-l_1.*l_2.*m_w.*r.^2.*cos(alpha-varphi_p).*2.0+l_1.*m_b.*r.^2.*y_G.*cos(alpha-varphi_p)+l_1.*m_b.*r.^2.*x_G.*sin(alpha-varphi_p).*2.0-l_2.*m_b.*r.^2.*y_G.*sin(alpha-varphi_p)).*(-1.0./2.0))./(l_1.*l_2.*r),(I_a.*l_2.^2.*sin(alpha-varphi_p).*2.0+I_b.*r.^2.*sin(alpha-varphi_p)+I_t.*r.^2.*sin(alpha-varphi_p).*2.0+l_1.^2.*m_b.*r.^2.*sin(alpha-varphi_p)+l_1.^2.*m_p.*r.^2.*sin(alpha-varphi_p)+l_2.^2.*m_w.*r.^2.*sin(alpha-varphi_p).*2.0+I_a.*l_1.*l_2.*cos(alpha-varphi_p).*2.0+m_b.*r.^2.*x_G.^2.*sin(alpha-varphi_p)+m_b.*r.^2.*y_G.^2.*sin(alpha-varphi_p)+l_1.*m_p.*r.^2.*y_F.*cos(alpha)+l_1.*m_p.*r.^2.*x_F.*sin(alpha)+l_1.*l_2.*m_b.*r.^2.*cos(alpha-varphi_p)+l_1.*l_2.*m_p.*r.^2.*cos(alpha-varphi_p)+l_1.*l_2.*m_w.*r.^2.*cos(alpha-varphi_p).*2.0+l_1.*m_b.*r.^2.*y_G.*cos(alpha-varphi_p)+l_1.*m_b.*r.^2.*x_G.*sin(alpha-varphi_p).*2.0+l_2.*m_b.*r.^2.*y_G.*sin(alpha-varphi_p))./(l_1.*l_2.*r.*2.0),-m_p.*(y_F.*cos(alpha)+x_F.*sin(alpha)),(I_a.*l_2.^2.*cos(alpha-varphi_p).*2.0+I_b.*r.^2.*cos(alpha-varphi_p)+I_t.*r.^2.*cos(alpha-varphi_p).*2.0+l_1.^2.*m_b.*r.^2.*cos(alpha-varphi_p)+l_1.^2.*m_p.*r.^2.*cos(alpha-varphi_p)+l_2.^2.*m_w.*r.^2.*cos(alpha-varphi_p).*2.0+m_b.*r.^2.*x_G.^2.*cos(alpha-varphi_p)+m_b.*r.^2.*y_G.^2.*cos(alpha-varphi_p)+I_a.*l_1.*l_2.*sin(alpha-varphi_p).*2.0+l_1.*m_p.*r.^2.*x_F.*cos(alpha)-l_1.*m_p.*r.^2.*y_F.*sin(alpha)+l_1.*m_b.*r.^2.*x_G.*cos(alpha-varphi_p).*2.0-l_2.*m_b.*r.^2.*y_G.*cos(alpha-varphi_p)+l_1.*l_2.*m_b.*r.^2.*sin(alpha-varphi_p)+l_1.*l_2.*m_p.*r.^2.*sin(alpha-varphi_p)+l_1.*l_2.*m_w.*r.^2.*sin(alpha-varphi_p).*2.0-l_1.*m_b.*r.^2.*y_G.*sin(alpha-varphi_p))./(l_1.*l_2.*r.*2.0),((I_a.*l_2.^2.*cos(alpha-varphi_p).*2.0+I_b.*r.^2.*cos(alpha-varphi_p)+I_t.*r.^2.*cos(alpha-varphi_p).*2.0+l_1.^2.*m_b.*r.^2.*cos(alpha-varphi_p)+l_1.^2.*m_p.*r.^2.*cos(alpha-varphi_p)+l_2.^2.*m_w.*r.^2.*cos(alpha-varphi_p).*2.0+m_b.*r.^2.*x_G.^2.*cos(alpha-varphi_p)+m_b.*r.^2.*y_G.^2.*cos(alpha-varphi_p)-I_a.*l_1.*l_2.*sin(alpha-varphi_p).*2.0+l_1.*m_p.*r.^2.*x_F.*cos(alpha)-l_1.*m_p.*r.^2.*y_F.*sin(alpha)+l_1.*m_b.*r.^2.*x_G.*cos(alpha-varphi_p).*2.0+l_2.*m_b.*r.^2.*y_G.*cos(alpha-varphi_p)-l_1.*l_2.*m_b.*r.^2.*sin(alpha-varphi_p)-l_1.*l_2.*m_p.*r.^2.*sin(alpha-varphi_p)-l_1.*l_2.*m_w.*r.^2.*sin(alpha-varphi_p).*2.0-l_1.*m_b.*r.^2.*y_G.*sin(alpha-varphi_p)).*(-1.0./2.0))./(l_1.*l_2.*r),m_p.*(x_F.*cos(alpha)-y_F.*sin(alpha)),(r.*(I_p+m_p.*x_F.^2+m_p.*y_F.^2+l_1.*m_p.*x_F.*cos(varphi_p)-l_2.*m_p.*y_F.*cos(varphi_p)-l_2.*m_p.*x_F.*sin(varphi_p)-l_1.*m_p.*y_F.*sin(varphi_p)))./(l_2.*2.0),(r.*(I_p+m_p.*x_F.^2+m_p.*y_F.^2+l_1.*m_p.*x_F.*cos(varphi_p)+l_2.*m_p.*y_F.*cos(varphi_p)+l_2.*m_p.*x_F.*sin(varphi_p)-l_1.*m_p.*y_F.*sin(varphi_p)).*(-1.0./2.0))./l_2,I_p+m_p.*x_F.^2+m_p.*y_F.^2],[3,3]);
C_bars2 = reshape([-(I_a.*(alpha_dot-varphi_dot_p).*(l_2.*cos(alpha-varphi_p)+l_1.*sin(alpha-varphi_p)))./(l_1.*r)+(r.*sin(alpha-varphi_p).*(alpha_dot-varphi_dot_p).*(l_1.*l_2.*m_w.*-2.0+l_2.*m_b.*x_G+l_1.*m_b.*y_G))./(l_1.*l_2.*2.0)-(r.*cos(alpha-varphi_p).*(alpha_dot-varphi_dot_p).*(I_b+I_t.*2.0+l_2.^2.*m_w.*2.0+m_b.*x_G.^2+m_b.*y_G.^2+l_1.*m_b.*x_G-l_2.*m_b.*y_G))./(l_1.*l_2.*2.0),(I_a.*(alpha_dot-varphi_dot_p).*(l_2.*cos(alpha-varphi_p)-l_1.*sin(alpha-varphi_p)))./(l_1.*r)-(r.*sin(alpha-varphi_p).*(alpha_dot-varphi_dot_p).*(l_1.*l_2.*m_w.*2.0-l_2.*m_b.*x_G+l_1.*m_b.*y_G))./(l_1.*l_2.*2.0)+(r.*cos(alpha-varphi_p).*(alpha_dot-varphi_dot_p).*(I_b+I_t.*2.0+l_2.^2.*m_w.*2.0+m_b.*x_G.^2+m_b.*y_G.^2+l_1.*m_b.*x_G+l_2.*m_b.*y_G))./(l_1.*l_2.*2.0),0.0,(I_a.*(alpha_dot-varphi_dot_p).*(l_1.*cos(alpha-varphi_p)-l_2.*sin(alpha-varphi_p)))./(l_1.*r)-(r.*cos(alpha-varphi_p).*(alpha_dot-varphi_dot_p).*(l_1.*l_2.*m_w.*-2.0+l_2.*m_b.*x_G+l_1.*m_b.*y_G))./(l_1.*l_2.*2.0)-(r.*sin(alpha-varphi_p).*(alpha_dot-varphi_dot_p).*(I_b+I_t.*2.0+l_2.^2.*m_w.*2.0+m_b.*x_G.^2+m_b.*y_G.^2+l_1.*m_b.*x_G-l_2.*m_b.*y_G))./(l_1.*l_2.*2.0),(I_a.*(alpha_dot-varphi_dot_p).*(l_1.*cos(alpha-varphi_p)+l_2.*sin(alpha-varphi_p)))./(l_1.*r)+(r.*cos(alpha-varphi_p).*(alpha_dot-varphi_dot_p).*(l_1.*l_2.*m_w.*2.0-l_2.*m_b.*x_G+l_1.*m_b.*y_G))./(l_1.*l_2.*2.0)+(r.*sin(alpha-varphi_p).*(alpha_dot-varphi_dot_p).*(I_b+I_t.*2.0+l_2.^2.*m_w.*2.0+m_b.*x_G.^2+m_b.*y_G.^2+l_1.*m_b.*x_G+l_2.*m_b.*y_G))./(l_1.*l_2.*2.0),0.0,(alpha_dot.*m_p.*r.*(l_2.*x_F.*cos(varphi_p)+l_1.*y_F.*cos(varphi_p)+l_1.*x_F.*sin(varphi_p)-l_2.*y_F.*sin(varphi_p)).*(-1.0./2.0))./l_2,(alpha_dot.*m_p.*r.*(-l_2.*x_F.*cos(varphi_p)+l_1.*y_F.*cos(varphi_p)+l_1.*x_F.*sin(varphi_p)+l_2.*y_F.*sin(varphi_p)))./(l_2.*2.0),0.0],[3,3]);
MIIKdots = reshape([-(l_2.*cos(alpha-varphi_p).*(alpha_dot-varphi_dot_p)+l_1.*sin(alpha-varphi_p).*(alpha_dot-varphi_dot_p))./(l_1.*r),(l_2.*cos(alpha-varphi_p).*(alpha_dot-varphi_dot_p)-l_1.*sin(alpha-varphi_p).*(alpha_dot-varphi_dot_p))./(l_1.*r),(cos(alpha-varphi_p).*(alpha_dot-varphi_dot_p))./l_1,(l_1.*cos(alpha-varphi_p).*(alpha_dot-varphi_dot_p)-l_2.*sin(alpha-varphi_p).*(alpha_dot-varphi_dot_p))./(l_1.*r),(l_1.*cos(alpha-varphi_p).*(alpha_dot-varphi_dot_p)+l_2.*sin(alpha-varphi_p).*(alpha_dot-varphi_dot_p))./(l_1.*r),(sin(alpha-varphi_p).*(alpha_dot-varphi_dot_p))./l_1,0.0,0.0,0.0],[3,3]);
Dmats = reshape([(l_2.*r.*cos(alpha-varphi_p)-l_1.*r.*sin(alpha-varphi_p))./(l_2.*2.0),(l_1.*r.*cos(alpha-varphi_p)+l_2.*r.*sin(alpha-varphi_p))./(l_2.*2.0),r./(l_2.*2.0),1.0,0.0,0.0,(l_2.*r.*cos(alpha-varphi_p)+l_1.*r.*sin(alpha-varphi_p))./(l_2.*2.0),((l_1.*r.*cos(alpha-varphi_p)-l_2.*r.*sin(alpha-varphi_p)).*(-1.0./2.0))./l_2,(r.*(-1.0./2.0))./l_2,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,1.0],[6,3]);

Msys2 = [M_bars2, zeros(3);
         -MIIKs, eye(3)];
     
%------------ Compute the equations ------------%

qdot = xs(7:12);
% (!)Warning: given that the u vector from the data is transposed we have to
% transpose it here again that is why in the forula it appears u' insted of
% directly u. We have to take into consideration that in this contex u is a
% 1x3 vector and we want a 3x1.

% Full equations including pivot force disturbances and fricction
% qdotdot = pinv(Msys2)*[u' + Dmats.'*tau_friction + Dmats.'*D_vec - C_bars2*qdot(1:3,1); MIIKdots*qdot(1:3,1)]; 

% Equations including only fricction
qdotdot = pinv(Msys2)*[u' + Dmats.'*tau_friction - C_bars2*qdot(1:3,1); MIIKdots*qdot(1:3,1)]; 

% Equations of motion without force disturbances or friction
% qdotdot = pinv(Msys2)*[u' - C_bars2*qdot(1:3,1); MIIKdots*qdot(1:3,1)];

xdot = [qdot; qdotdot];

%------------ Matrices to compute IMU readings ------------%

% Now we configure the necesary matrices 
% Rzmat22 = reshape([cos(alpha),sin(alpha),-sin(alpha),cos(alpha)],[2,2]);
IRzmat22 = reshape([cos(alpha)./(cos(alpha).^2+sin(alpha).^2),-sin(alpha)./(cos(alpha).^2+sin(alpha).^2),sin(alpha)./(cos(alpha).^2+sin(alpha).^2),cos(alpha)./(cos(alpha).^2+sin(alpha).^2)],[2,2]);

% This matrix is the time derivative of the inverse rotation matrix (only needed if
% we use the old way of computing IMU readings)
% dIRzmat22 = reshape([-alpha_dot.*sin(alpha),-alpha_dot.*cos(alpha),alpha_dot.*cos(alpha),-alpha_dot.*sin(alpha)],[2,2]);

%-------------------- Output equations --------------------%

% Output equation with optional bias (scaling factor)
% yout = [xs(9) + Abias*xs(9);       % alpha_dot
%         qdotdot(1) + Xbias*qdotdot(1);  % x_dotdot
%         qdotdot(2) + Ybias*qdotdot(2)]; % y_dotdot

% Old equations of IMU readings (using time derivatives of the inverse
% rotation matrix)
% xyIMU = dIRzmat22*[x_dot;y_dot] + IRzmat22*[qdotdot(1);qdotdot(2)];

% New equations of IMU readings (projection of the acceleration vector to another base)
xyIMU = IRzmat22*[qdotdot(1);qdotdot(2)];

% Output equation now with optional bias and IMU readings
% yout = [xs(9) + Abias*xs(9);        % alpha_dot
%         xyIMU(1) + Xbias*xyIMU(1);  % x_dotdot
%         xyIMU(2) + Ybias*xyIMU(2)]; % y_dotdot 

% Output equation using only IMU readings
yout = [xs(9);     % alpha_dot
        xyIMU(1);  % x_dotdot
        xyIMU(2)]; % y_dotdot 

end


