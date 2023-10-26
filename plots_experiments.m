%% Leader Profiles

clearvars;

results_folder = "final_plots";

%%Plot leader velocity
load("results/real_traj/init3/3av_5.6667hv_1/leader.mat", "leader")
load("results/real_traj/init3/3av_5.6667hv_1/auxadata_10.mat")
vl = leader.v(auxdata.time);
figure()
plot(auxdata.time, vl, "color", '#5E5E5E', "linewidth", 1)
title("Leader's Trajectory - Velocity", "fontsize", 12)
xlabel("Time (s)")
ylabel("Velocity (m/s)")
savefig(results_folder + '/leader/velocity_leader.fig')
close()

%%Plot leader acceleration  
al = diff(vl)./diff(auxdata.time);
al = [al; al(end)];
figure()
plot(auxdata.time, al, "color", '#5E5E5E', "linewidth", 1)
title("Leader's Trajectory - Acceleration", "fontsize", 12)
xlabel("Time (s)")
ylabel("Acceleration (m/s^2)")
savefig(results_folder + '/leader/accel_leader.fig')
close()

%% Plot Optimized AV Behavior - Penetration Experiments 
clearvars;

% Set of colors
colors = ["#e6194B", "#f58231", "#ffe119", "#bfef45", "#3cb44b"];

% UPDATE PATH
load("results/real_traj/init3/0av/leader.mat", "leader")
load("results/real_traj/init3/0av/auxadata_10.mat", "auxdata")
load("results/real_traj/init3/0av/U_10.mat", "U_star")

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
fig_filename = results_folder + "/velocity";
savefig(fig_filename + ".fig")
fig=openfig(fig_filename + ".fig",'new','invisible');
saveas(fig,fig_filename + ".png",'png');
close(fig);


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
fig_filename = results_folder + "/acceleration";
savefig(fig_filename + ".fig")
fig=openfig(fig_filename + ".fig",'new','invisible');
saveas(fig,fig_filename + ".png",'png');
close(fig);

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
ylim([-5 125]);
title("Penetration Experiment - " + length(auxdata.Ia) + "AVs", "fontsize", 12)
xlabel("Time (s)", "fontsize", 12)
ylabel("Headway (m)", "fontsize", 12)
legend(plots_list, legend_arr)
fig_filename = results_folder + "/headway";
savefig(fig_filename + ".fig")
fig=openfig(fig_filename + ".fig",'new','invisible');
saveas(fig,fig_filename + ".png",'png');
close(fig);

%% Plot Greedy Experiments

clearvars;

% Set of colors
colors = ["#e6194B", "#f58231", "#ffe119", "#bfef45", "#3cb44b"];

% UPDATE PATH
platoon_name = "5av_2.4hv_1"
path = "results/real_traj/greedy/new_greedy/" + platoon_name + "/";
load(path+"leader.mat", "leader")
load(path+"auxadata_10", "auxdata")
if platoon_name == "1av_19hv"
    load(path+"U_10.mat", "U_1av")
    U_star = U_1av;
else
    load(path+"U_10.mat", "U_star")
end

%Path to store results
results_folder = "final_plots/new_greedy_experiments/" + platoon_name;
if not(isfolder(results_folder))
    mkdir(results_folder)
end

% Get X_star, V_star, A_star
[X_star, V_star, A_star] = system_solve(U_star, auxdata, leader);
display(trapz(auxdata.time, sum(A_star.^2, 2)))
display(sum(trapz(auxdata.time, simplified_fuel_model(V_star, A_star, 'RAV4'))))

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
title("Greedy Experiment - " + length(auxdata.Ia) + "AVs", "fontsize", 12)
xlabel("Time (s)")
ylabel("Velocity (m/s)")
legend(plots_list, legend_arr)
fig_filename = results_folder + "/velocity";
savefig(fig_filename + ".fig")
fig=openfig(fig_filename + ".fig",'new','invisible');
saveas(fig,fig_filename + ".png",'png');
close(fig);


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
title("Greedy Experiment - " + length(auxdata.Ia) + "AVs", "fontsize", 12)
xlabel("Time (s)")
ylabel("Acceleration (m/s^2)")
legend(plots_list, legend_arr)
fig_filename = results_folder + "/acceleration";
savefig(fig_filename + ".fig")
fig=openfig(fig_filename + ".fig",'new','invisible');
saveas(fig,fig_filename + ".png",'png');
close(fig);

%Plot Optimal Solution Headways - Greedy Experiments
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
ylim([-5 125]);
title("Greedy Experiment - " + length(auxdata.Ia) + "AVs", "fontsize", 12)
xlabel("Time (s)", "fontsize", 12)
ylabel("Headway (m)", "fontsize", 12)
legend(plots_list, legend_arr)
fig_filename = results_folder + "/headway";
savefig(fig_filename + ".fig")
fig=openfig(fig_filename + ".fig",'new','invisible');
saveas(fig,fig_filename + ".png",'png');
close(fig);

%% Plot 0AV 
clearvars;

% Set of colors
colors = ["#e6194B", "#f58231", "#ffe119", "#bfef45", "#3cb44b"];

% UPDATE PATH
load("results/real_traj/init3/1av_19hv_1/leader.mat", "leader")
load("results/real_traj/init3/1av_19hv_1/auxadata_10.mat", "auxdata")

platoon_name = "0av_20hv";
results_in = "results/real_traj/init3" + "/" + platoon_name + "_1/"; 

%Path to store results
results_folder = "final_plots/penetration_experiments/" + platoon_name;
if not(isfolder(results_folder))
    mkdir(results_folder)
end

% Get X_star, V_star, A_star
load("results/real_traj/init3/0av/A_0av.mat")
load("results/real_traj/init3/0av/V_0av.mat")
load("results/real_traj/init3/0av/X_0av.mat")

auxdata.platoon = [zeros(1,20)];

%%Plot leader leader velocity and optimized AV velocities  
vl = leader.v(auxdata.time);
legend_arr = ["Leader"];
figure()
leader_p = plot(auxdata.time, vl, "color", 'black', "linewidth", 4);
plots_list = [leader_p];
color_iter = 1;
hold on 
p = 0;
for i=1:auxdata.len_platoon
    if auxdata.platoon(i) == 0
        p = plot(auxdata.time, V1av(:, i), "linewidth", 0.5, "Color","blue");
        hold on
    else
        hold on
    end
end
plots_list(end+1) = p;
legend_arr(end+1) = "HVs";
title("Penetration Experiment - " + length(auxdata.Ia) + "AVs", "fontsize", 12)
xlabel("Time (s)")
ylabel("Velocity (m/s)")
legend(plots_list, legend_arr)
fig_filename = results_folder + "/velocity";
savefig(fig_filename + ".fig")
fig=openfig(fig_filename + ".fig",'new','invisible');
saveas(fig,fig_filename + ".png",'png');
close(fig);


%%Plot leader acceleration with optimized AV Acceleration
al = diff(vl)./diff(auxdata.time);
al = [al; al(end)];
legend_arr = ["Leader"];
figure()
leader_p = plot(auxdata.time, al, "color", 'black', "linewidth", 4);
plots_list = [leader_p];
color_iter = 1;
hold on 
p = 0;
for i=1:auxdata.len_platoon
    if auxdata.platoon(i) == 0
        p = plot(auxdata.time, A1av(:, i), "linewidth", 0.5, "Color","blue");
        hold on
    else
        hold on
    end
end
plots_list(end+1) = p;
legend_arr(end+1) = "HVs";
title("Penetration Experiment - " + length(auxdata.Ia) + "AVs", "fontsize", 12)
xlabel("Time (s)")
ylabel("Acceleration (m/s^2)")
legend(plots_list, legend_arr)
fig_filename = results_folder + "/acceleration";
savefig(fig_filename + ".fig")
fig=openfig(fig_filename + ".fig",'new','invisible');
saveas(fig,fig_filename + ".png",'png');
close(fig);
