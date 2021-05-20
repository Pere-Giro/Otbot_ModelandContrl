Geom_Vec=from_file('geom_data_otbot.txt');

m.I_b=Geom_Vec(1,1);
m.I_p=Geom_Vec(2,1);
m.I_a=Geom_Vec(3,1);
m.I_t=Geom_Vec(4,1);
m.l_1=Geom_Vec(5,1);
m.l_2=Geom_Vec(6,1);
m.m_b=Geom_Vec(7,1);
m.m_w=Geom_Vec(8,1);
m.m_p=Geom_Vec(9,1);
m.x_G=Geom_Vec(10,1);
m.y_G=Geom_Vec(11,1);
m.x_F=Geom_Vec(12,1);
m.y_F=Geom_Vec(13,1);
m.r=Geom_Vec(14,1);

q.x = -0.3;
q.y = 1.5;
q.alpha = pi/2;
q.varphi_r = 0;
q.varphi_l = 0;
q.varphi_p = 0;

draw_otbot_editing(m,q)

%%

rectangle('Position',[0.5 0.5 1 1])

ABC = [0, 1, 1, 0;
       0, 0, 1, 1;
       0, 0, 0, 0];
   
multi_circles(ABC,0.1/6,'r')

axis equal