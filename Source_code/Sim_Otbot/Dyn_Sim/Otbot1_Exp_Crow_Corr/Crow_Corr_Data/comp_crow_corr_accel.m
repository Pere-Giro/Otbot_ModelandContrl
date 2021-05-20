%% Loading the data

load("crow_corr_ts")

u_actions = [];

for t = 0:0.001:9.99
    u_f = getdatasamples(resample(uactionts,t,'linear'),1);
    u_actions = cat(2,u_actions,u_f);
end

%% Compute the acceleration vector
load("sm_struc")
pddotvec = [];

for i = 1:9991
    alpha = sstates(3,i);
    varphi_p = sstates(7,i);
    x_dot = sstates(4,i);
    y_dot = sstates(5,i);
    alpha_dot = sstates(6,i);
    
    MIIKs = sm.MIIKmatrix(alpha,varphi_p);
    
    phi_dot = MIIKs*[x_dot;y_dot;alpha_dot];
    
    varphi_dot_p = phi_dot(3,1);
    
    M_bars2 = sm.M_barmatrix2(alpha,varphi_p);
    C_bars2 = sm.C_barmatrix2(alpha,alpha_dot,varphi_p,varphi_dot_p);
    MIIKdots = sm.MIIKdotmatrix(alpha,alpha_dot,varphi_p,varphi_dot_p);
    Dmats = sm.Deltamatrix(alpha,varphi_p);

    Msys2 = [M_bars2, zeros(3);
             -MIIKs, eye(3)];
         
    qdot = [x_dot;y_dot;alpha_dot;phi_dot];
    
    u_f = u_actions(:,i);
    qdotdot = pinv(Msys2)*[u_f - C_bars2*qdot(1:3,1); MIIKdots*qdot(1:3,1)];
    
    pddot = qdotdot(1:3,1);
    
    pddotvec = cat(2,pddotvec,pddot);
end

%% Create the timeseries object

pddotts = timeseries(pddotvec,[0:0.001:9.99]);

save("pddot_ts.mat","pddotts")

    
    
    
    
    
    