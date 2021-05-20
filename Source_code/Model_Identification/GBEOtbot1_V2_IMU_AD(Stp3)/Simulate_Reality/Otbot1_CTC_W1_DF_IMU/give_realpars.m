function real_pars = give_realpars(m,Biasalphadot_Flag, Biasxddot_Flag, Biasyddot_Flag, IMU)
%GIVE_REALPARS Summary of this function goes here
%   This function returs the real pars vector depending on the flags
%   settings

Abias = IMU.Abias;
Xbias = IMU.Xbias;
Ybias = IMU.Ybias;

if strcmp(Biasalphadot_Flag,'NO') && strcmp(Biasxddot_Flag,'NO') && strcmp(Biasyddot_Flag,'NO')
    real_pars = [realpfunc(m);
                 0;
                 0;
                 0];
    
elseif strcmp(Biasalphadot_Flag,'NO') && strcmp(Biasxddot_Flag,'NO') && strcmp(Biasyddot_Flag,'YES')
    real_pars = [realpfunc(m);
                 0;
                 0;
                 Ybias];
    
elseif strcmp(Biasalphadot_Flag,'NO') && strcmp(Biasxddot_Flag,'YES') && strcmp(Biasyddot_Flag,'NO')
    real_pars = [realpfunc(m);
                 0;
                 Xbias;
                 0];
    
elseif strcmp(Biasalphadot_Flag,'NO') && strcmp(Biasxddot_Flag,'YES') && strcmp(Biasyddot_Flag,'YES')
    real_pars = [realpfunc(m);
                 0;
                 Xbias;
                 Ybias];
    
elseif strcmp(Biasalphadot_Flag,'YES') && strcmp(Biasxddot_Flag,'NO') && strcmp(Biasyddot_Flag,'NO')
    real_pars = [realpfunc(m);
                 Abias;
                 0;
                 0];

elseif strcmp(Biasalphadot_Flag,'YES') && strcmp(Biasxddot_Flag,'NO') && strcmp(Biasyddot_Flag,'YES')
    real_pars = [realpfunc(m);
                 Abias;
                 0;
                 Ybias];
    
elseif strcmp(Biasalphadot_Flag,'YES') && strcmp(Biasxddot_Flag,'YES') && strcmp(Biasyddot_Flag,'NO')
    real_pars = [realpfunc(m);
                 Abias;
                 Xbias;
                 0];
    
elseif strcmp(Biasalphadot_Flag,'YES') && strcmp(Biasxddot_Flag,'YES') && strcmp(Biasyddot_Flag,'YES')
    real_pars = [realpfunc(m);
                 Abias;
                 Xbias;
                 Ybias];
    
else
    disp('WARNING: This Biasalphadot_Flag, Biasxddot_Flag and Biasyddot_Flag combination does not exist, real pars vector builded considering no scaling factor error')
    real_pars = [realpfunc(m);
                 0;
                 0;
                 0];
    
end



end

