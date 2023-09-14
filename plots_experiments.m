%% Leader Profiles

results_folder = "final_plots";

%%Plot leader velocity
load("results/real_traj/init3/3av_5.6667hv_1/leader.mat", "leader")
vl = leader.v(auxdata.time);
figure()
plot(auxdata.time, vl, "color", '#5E5E5E', "linewidth", 1)
title("Leader's Trajectory - Velocity", "fontsize", 12)
xlabel("Time (s)")
ylabel("Velocity (m/s)")
savefig(results_folder + '/leader/velocity_leader.fig')

%%Plot leader acceleration  
al = diff(vl)./diff(auxdata.time);
al = [al; al(end)];
figure()
plot(auxdata.time, al, "color", '#5E5E5E', "linewidth", 1)
title("Leader's Trajectory - Acceleration", "fontsize", 12)
xlabel("Time (s)")
ylabel("Acceleration (m/s^2)")
savefig(results_folder + '/leader/accel_leader.fig')

%% Plot Optimized AV Behavior - Penetration Experiments 
clearvars;

% Set of colors
colors = ["#e6194B", "#f58231", "#ffe119", "#bfef45", "#3cb44b"];

% UPDATE PATH
load("results/real_traj/init3/2av_9hv_1/leader.mat", "leader")
load("results/real_traj/init3/2av_9hv_1/auxadata_10.mat", "auxdata")
load("results/real_traj/init3/2av_9hv_1/U_10.mat", "U_star")

platoon_name = length(auxdata.Ia) + "av_" + length(auxdata.Ih)/length(auxdata.Ia) + "hv";
results_in = "results/real_traj/init3" + "/" + platoon_name + "_1/"; 

%Path to store results
results_folder = "final_plots/penetration_experiments/" + platoon_name;
if not(isfolder(results_folder))
    mkdir(results_folder)
end

% Get X_star, V_star, A_star
[X_star, V_star, A_star] = system_solve(U_star, auxdata, leader); 

%%Plot leader leader velocity and optimized AV velocities  
vl = leader.v(auxdata.time);
legend_arr = ["Leader"];
figure()
leader_p = plot(auxdata.time, vl, "color", 'black', "linewidth", 1);
plots_list = [leader_p];
color_iter = 1;
hold on 
p = 0;
for i=1:auxdata.len_platoon
    if auxdata.platoon(i) == 0
        p = plot(auxdata.time, V_star(:, i), "linewidth", 0.5, "Color","blue");
        hold on
    else
        p = plot(auxdata.time, V_star(:, i), "linewidth", 4, "Color", colors(color_iter));
        color_iter = color_iter+1;
        plots_list(end+1) = p;
        legend_arr(end+1) = "V" + string(i) + " (AV)";
        hold on
    end
end
plots_list(end+1) = p;
legend_arr(end+1) = "HVs";
title("Penetration Experiment - " + length(auxdata.Ia) + "AVs", "fontsize", 12)
xlabel("Time (s)")
ylabel("Velocity (m/s)")
legend(plots_list, legend_arr)
savefig(results_folder + "/velocity.fig")


%%Plot leader acceleration with optimized AV Acceleration
al = diff(vl)./diff(auxdata.time);
al = [al; al(end)];
legend_arr = ["Leader"];
figure()
leader_p = plot(auxdata.time, al, "color", 'black', "linewidth", 1);
plots_list = [leader_p];
color_iter = 1;
hold on 
p = 0;
for i=1:auxdata.len_platoon
    if auxdata.platoon(i) == 0
        p = plot(auxdata.time, A_star(:, i), "linewidth", 0.5, "Color","blue");
        hold on
    else
        p = plot(auxdata.time, A_star(:, i), "linewidth", 4, "Color", colors(color_iter));
        color_iter = color_iter+1;
        plots_list(end+1) = p;
        legend_arr(end+1) = "V" + string(i) + " (AV)";
        hold on
    end
end
plots_list(end+1) = p;
legend_arr(end+1) = "HVs";
title("Penetration Experiment - " + length(auxdata.Ia) + "AVs", "fontsize", 12)
xlabel("Time (s)")
ylabel("Acceleration (m/s^2)")
legend(plots_list, legend_arr)
savefig(results_folder + "/acceleration.fig")

%Plot Optimal Solution Headways - Penetration Experiments
% Set of colors
colors = ["#e6194B", "#f58231", "#ffe119", "#bfef45", "#3cb44b"];
legend_arr = ["V1 (AV) -> Leader"];
Xl = [leader.x(auxdata.time), X_star]; 
headway = Xl(:,auxdata.Ia) - X_star(:, auxdata.Ia) - auxdata.l;
figure()
leader_p = plot(auxdata.time, headway(:, 1), "color", "#e6194B", "linewidth", 2);
plots_list = [leader_p];
color_iter = 2;
hold on 
p = 0;
for i=2:size(headway, 2)
    p = plot(auxdata.time, headway(:, i), "linewidth", 2, "Color", colors(color_iter));
    color_iter = color_iter+1;
    plots_list(end+1) = p;
    legend_arr(end+1) = "V" + string(4*i - 3) + " (AV)" + " -> " + "V" + string(4*i-3-1) + " (HV)";
    hold on
end
plot(auxdata.time, yline(auxdata.d_min), "color", 'black', "linewidth", 5, "LineStyle", "--")
hold on 
plot(auxdata.time, yline(auxdata.d_max), "color", 'black', "linewidth", 5, "LineStyle", "--")
hold on
title("Penetration Experiment - " + length(auxdata.Ia) + "AVs", "fontsize", 12)
xlabel("Time (s)", "fontsize", 12)
ylabel("Headway (m)", "fontsize", 12)
legend(plots_list, legend_arr)
savefig(results_folder + "/headway.fig")

%% Plot Greedy Experiments

clearvars;

% Set of colors
colors = ["#e6194B", "#f58231", "#ffe119", "#bfef45", "#3cb44b"];

% UPDATE PATH
path = "results/real_traj/greedy/init1/5av_3hv_1/";
load(path+"leader.mat", "leader")
load(path+"auxadata_10.mat", "auxdata")
load(path+"U_10.mat", "U_star")

platoon_name = length(auxdata.Ia) + "av_" + length(auxdata.Ih)/length(auxdata.Ia) + "hv";

%Path to store results
results_folder = "final_plots/greedy_experiments/" + platoon_name;
if not(isfolder(results_folder))
    mkdir(results_folder)
end

% Get X_star, V_star, A_star
[X_star, V_star, A_star] = system_solve(U_star, auxdata, leader); 

%%Plot leader leader velocity and optimized AV velocities  
vl = leader.v(auxdata.time);
legend_arr = ["Leader"];
figure()
leader_p = plot(auxdata.time, vl, "color", 'black', "linewidth", 0.5);
plots_list = [leader_p];
color_iter = 1;
hold on 
p = 0;
for i=1:auxdata.len_platoon
    if auxdata.platoon(i) == 0
        p = plot(auxdata.time, V_star(:, i), "linewidth", 1, "LineStyle", ":", "Color","blue");
        hold on
    else
        p = plot(auxdata.time, V_star(:, i), "linewidth", 4, "Color", colors(color_iter));
        color_iter = color_iter+1;
        plots_list(end+1) = p;
        legend_arr(end+1) = "V" + string(i) + " (AV)";
        hold on
    end
end
plots_list(end+1) = p;
legend_arr(end+1) = "HVs";
title("Penetration Experiment - " + length(auxdata.Ia) + "AVs", "fontsize", 12)
xlabel("Time (s)")
ylabel("Velocity (m/s)")
legend(plots_list, legend_arr)
savefig(results_folder + "/velocity.fig")


%%Plot leader acceleration with optimized AV Acceleration
al = diff(vl)./diff(auxdata.time);
al = [al; al(end)];
legend_arr = ["Leader"];
figure()
leader_p = plot(auxdata.time, al, "color", 'black', "linewidth", 0.5);
plots_list = [leader_p];
color_iter = 1;
hold on 
p = 0;
for i=1:auxdata.len_platoon
    if auxdata.platoon(i) == 0
        p = plot(auxdata.time, A_star(:, i), "linewidth", 0.75, "LineStyle", ":", "Color","blue");
        hold on
    else
        p = plot(auxdata.time, A_star(:, i), "linewidth", 4, "Color", colors(color_iter));
        color_iter = color_iter+1;
        plots_list(end+1) = p;
        legend_arr(end+1) = "V" + string(i) + " (AV)";
        hold on
    end
end
plots_list(end+1) = p;
legend_arr(end+1) = "HVs";
title("Penetration Experiment - " + length(auxdata.Ia) + "AVs", "fontsize", 12)
xlabel("Time (s)")
ylabel("Acceleration (m/s^2)")
legend(plots_list, legend_arr)
savefig(results_folder + "/acceleration.fig")

%Plot Optimal Solution Headways - Penetration Experiments
% Set of colors
colors = ["#e6194B", "#f58231", "#ffe119", "#bfef45", "#3cb44b"];
legend_arr = ["V1 -> Leader"];
Xl = [leader.x(auxdata.time), X_star]; 
headway = Xl(:,auxdata.Ia) - X_star(:, auxdata.Ia) - auxdata.l;
figure()
leader_p = plot(auxdata.time, headway(:, 1), "color", "#e6194B", "linewidth", 3);
plots_list = [leader_p];
color_iter = 2;
hold on 
p = 0;
for i=2:size(headway, 2)
    p = plot(auxdata.time, headway(:, i), "linewidth", 3, "Color", colors(color_iter));
    color_iter = color_iter+1;
    plots_list(end+1) = p;
    legend_arr(end+1) = "V" + string(4*i - 3) + " -> " + "V" + string(4*i-3-1);
    hold on
end
plot(auxdata.time, yline(auxdata.d_min), "color", 'black', "linewidth", 5, "LineStyle", "--")
hold on 
plot(auxdata.time, yline(auxdata.d_max), "color", 'black', "linewidth", 5, "LineStyle", "--")
hold on
title("Penetration Experiment - " + length(auxdata.Ia) + "AVs", "fontsize", 12)
xlabel("Time (s)", "fontsize", 12)
ylabel("Headway (m)", "fontsize", 12)
legend(plots_list, legend_arr)
savefig(results_folder + "/headway.fig")