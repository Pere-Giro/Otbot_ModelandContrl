%% First load the required data

load("crow_corr_ts.mat")

h = 0.001;

times_ts = 0:h:10;

x_series = sstates(1,:);
y_series = sstates(2,:);
alpha_series = sstates(3,:);

x_dot_series = diff(x_series)/h;
y_dot_series = diff(y_series)/h;
alpha_dot_series = diff(alpha_series)/h;

x_ddot_series = diff(x_dot_series)/h;
y_ddot_series = diff(y_dot_series)/h;
alpha_ddot_series = diff(alpha_dot_series)/h;

%% Create the timeseries

CTC_vec = [x_series(1:9999);
           y_series(1:9999);
           alpha_series(1:9999);
           x_dot_series(1:9999);
           y_dot_series(1:9999);
           alpha_dot_series(1:9999);
           x_ddot_series;
           y_ddot_series;
           alpha_ddot_series];

CTCd_ts = timeseries(CTC_vec,times_ts(1:9999));
       
% x_ts = timeseries(x_series,times_ts);
% y_ts = timeseries(y_series,times_ts);
% alpha_ts = timeseries(alpha_series,times_ts);
% 
% dx_ts = timeseries(y_dot_series,times_ts(1:10000));
% dy_ts = timeseries(y_dot_series,times_ts(1:10000));
% dalpha_ts = timeseries(y_dot_series,times_ts(1:10000));
% 
% ddx_ts = timeseries(y_ddot_series,times_ts(1:9999));
% ddy_ts = timeseries(y_ddot_series,times_ts(1:9999));
% ddalpha_ts = timeseries(y_ddot_series,times_ts(1:9999));

%% Saving the series
% save("crow_corr_diff_traj.mat","x_ts","y_ts","alpha_ts","dx_ts","dy_ts","dalpha_ts","ddx_ts","ddy_ts","ddalpha_ts")
save("crow_corr_diff_traj.mat","CTCd_ts")

