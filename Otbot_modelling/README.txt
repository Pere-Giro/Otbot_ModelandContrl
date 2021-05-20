

Contents of the folders:

------------------------
Initial_pars
------------------------

These are the mlx files used to create the kinematic and dynamic models as explained in Maxime's thesis. 

Here, the chassis is separeated into three bodies to compute the translational kinetic energy: the left wheel, the right wheel, and the chassis main body.

Thus we here use:

	m_b = Mass of the main body of the chassis
	I_b = Moment of inertia of the chassis main body relative to its c.o.m. G
	(x_G,y_G) = Coords. of the c.o.m. of the chassis main body
	m_w = Mass of a wheel
	I_w = Vertical moment of inertia of a wheel at its c.o.m


------------------------
New_pars
------------------------

These are the mlx files used to create the kinematic and dynamic models as explained in Pere's thesis. 

Here the chassis main body and the wheels are viewed as a single body when computing the translational kinetic energy. This single body is called the "chassis".

Thus we here use:

	m_c = Mass of the chassis
	I_c = Moment of inertia of the chassis
	(x_B,y_B) = Coords. of the c.o.m. of the chassis

So this involves two less parameters in comparison to the initial set of parameters, which is more convenient from the point of view of system identification. Note that it is not possible to separately identify m_b and m_w by moving the robot in the plane. We can only identify the overall chassis mass m_c. A similar thing occurs for the moments of inertia I_b and I_w. Only I_c can be identified.



