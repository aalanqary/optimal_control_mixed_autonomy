ss = dir("data_v2_preprocessed_west");
i = 48;
ss(i).name
data = readtable("data_v2_preprocessed_west/"+ss(i).name);
time = data.Time;
time = time - time(1);
T = floor(time(end));
auxdata.dt = 0.1;
auxdata.udt = 1; 
auxdata.utime = (0:auxdata.udt:T)';
auxdata.time = (0:auxdata.dt:T)';

vl = data.Velocity * (1000/3600); 
plot(vl)