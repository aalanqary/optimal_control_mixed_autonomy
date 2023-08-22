%% Leader Profiles (FROM ARWA)

%%Plot leader velocity  
vl = auxdata.vl(auxdata.time);
figure()
plot(auxdata.time, vl, "color", '#5E5E5E', "linewidth", 1)
title("Leader's trajectory", "fontsize", 12)
xlabel("Time (s)")
ylabel("Velocity (m/s)")

%%Plot leader acceleration  
vl = auxdata.vl(auxdata.time);
al = diff(vl)./diff(auxdata.time);
al = [al; al(end)];
figure()
plot(auxdata.time, al, "color", '#5E5E5E', "linewidth", 1)
title("Leader's trajectory", "fontsize", 12)
xlabel("Time (s)")
ylabel("Acceleration (m/s^2)")

%% Initial AV Profiles
%%Plot leader position with first AV position -> Show initialization  

legend_arr = ["Leader"];
for k=1:auxdata.len_platoon
    legend_arr(end+1) = "V" + string(k);
end

xl = auxdata.xl(auxdata.time);
figure()
plot(auxdata.time, xl, "color", '#5E5E5E', "linewidth", 1)
hold on 
for i=1:auxdata.len_platoon
    plot(auxdata.time, X0(:, i), "linewidth", 2, "LineStyle", "--")
    hold on
end
title("Platoon Initialization [" + string(strjoin(string(auxdata.platoon))) + "] Simulation of Position", "fontsize", 12)
xlabel("Time (s)")
ylabel("Position (m)")
legend(legend_arr)

%%Plot leader velocity with first AV velocity -> Show initialization  
vl = auxdata.vl(auxdata.time);
figure()
plot(auxdata.time, vl, "color", '#5E5E5E', "linewidth", 1)
hold on 
for i=1:auxdata.len_platoon
    plot(auxdata.time, V0(:, i), "linewidth", 2, "LineStyle", "--")
    hold on
end
title("Platoon Initialization [" + string(strjoin(string(auxdata.platoon))) + "] Simulation of Velocity", "fontsize", 12)
xlabel("Time (s)")
ylabel("Velocity (m/s)")
legend(legend_arr)
%%Plot leader acceleration with first AV acceleration
vl = auxdata.vl(auxdata.time);
al = diff(vl)./diff(auxdata.time);
al = [al; al(end)];
figure()
plot(auxdata.time, al, "color", '#5E5E5E', "linewidth", 1)
hold on 
for i=1:auxdata.len_platoon
    plot(auxdata.time, A0(:, i), "linewidth", 2, "LineStyle", "--")
    hold on
end
title("Platoon Initialization [" + string(strjoin(string(auxdata.platoon))) + "] Simulation of Acceleration", "fontsize", 12)
xlabel("Time (s)")
ylabel("Acceleration (m/s^2)")
legend(legend_arr)

%% Plot headway feasibility of initial solution
% Plot min headway feasibility
figure()
xl = auxdata.xl(auxdata.time);
plot(auxdata.time, xl - X0(:, 1) - auxdata.l, "color", 'red', "linewidth", 1)

% hold on 
% plot(auxdata.time, auxdata.d_min * V0, "color", 'black', "linewidth", 0.5, "LineStyle", "--")
% plot(auxdata.time, xl - X0(:, 1) - auxdata.l, "color", '#EE220C', "linewidth", 1)
title("First AV -> Leader Headway", "fontsize", 12)
xlabel("Time (s)", "fontsize", 12)
ylabel("Headway (m)", "fontsize", 12)

% Ensure min headway feasibility

% Plot max headway feasibility

%% Calculate Initial System Energy Expenditure
% Calculate System Acceleration
init_system_accel = sum(trapz(auxdata.time, A0.^2))
% Calculate System Energy Expense
init_system_energy = sum(trapz(auxdata.time, simplified_fuel_model(V0, A0, 'RAV4')))

%% Plot Optimized AV Behavior

legend_arr = ["Leader"];
for k=1:auxdata.len_platoon
    if auxdata.platoon(k) == 1
        legend_arr(end+1) = "V" + string(k) + " (AV)";
    else
        legend_arr(end+1) = "V" + string(k);
    end
end

%%Plot leader position with AV position 
xl = auxdata.xl(auxdata.time);
figure()
plot(auxdata.time, xl, "color", '#5E5E5E', "linewidth", 1)
hold on 
for i=1:auxdata.len_platoon
    if auxdata.platoon(i) == 0
        plot(auxdata.time, X_star(:, i), "linewidth", 2, "LineStyle", "--")
        hold on
    else
        plot(auxdata.time, X_star(:, i), "linewidth", 2)
        hold on
    end
end
title("Mixed Autonomy Platoon [" + string(strjoin(string(auxdata.platoon))) + "] Simulation of Position", "fontsize", 12)
xlabel("Time (s)")
ylabel("Position (m)")
legend(legend_arr)

%%Plot leader leader velocity and optimized AV velocities  
vl = auxdata.vl(auxdata.time);
figure()
plot(auxdata.time, vl, "color", '#5E5E5E', "linewidth", 1)
hold on 
for i=1:auxdata.len_platoon
    if auxdata.platoon(i) == 0
        plot(auxdata.time, V_star(:, i), "linewidth", 2, "LineStyle", "--")
        hold on
    else
        plot(auxdata.time, V_star(:, i), "linewidth", 2)
        hold on
    end
end
title("Mixed Autonomy Platoon [" + string(strjoin(string(auxdata.platoon))) + "] Simulation of Velocity", "fontsize", 12)
xlabel("Time (s)")
ylabel("Velocity (m/s)")
legend(legend_arr)


%%Plot leader acceleration with optimized AV Acceleration
vl = auxdata.vl(auxdata.time);
al = diff(vl)./diff(auxdata.time);
al = [al; al(end)];
figure()
plot(auxdata.time, al, "color", '#5E5E5E', "linewidth", 1)
hold on 
for i=1:auxdata.len_platoon
    if auxdata.platoon(i) == 0
        plot(auxdata.time, A_star(:, i), "linewidth", 2, "LineStyle", "--")
        hold on
    else
        plot(auxdata.time, A_star(:, i), "linewidth", 2)
        hold on
    end
end
title("Mixed Autonomy Platoon [" + string(strjoin(string(auxdata.platoon))) + "] Simulation of Acceleration", "fontsize", 12)
xlabel("Time (s)")
ylabel("Acceleration (m/s^2)")
legend(legend_arr)

%% Plot Optimal Solution Headways
% Plot Optimal Solution Headway for AV1
legend_hw = ["V1 -> Leader"];
for k=1:auxdata.len_platoon-1
    legend_hw(end+1) = "V" + string(k+1) + " -> " + "V" + string(k);
end
figure()
xl = auxdata.xl(auxdata.time);
plot(auxdata.time, xl - X_star(:, 1) - auxdata.l, "color", '#EE220C', "linewidth", 1)
hold on 
for i=1:auxdata.len_platoon-1
    plot(auxdata.time, X_star(:, i) - X_star(:, i+1), "linewidth", 2)
    hold on
end
plot(auxdata.time, yline(auxdata.d_min), "color", 'black', "linewidth", 3, "LineStyle", "--")
hold on 
plot(auxdata.time, yline(auxdata.d_max), "color", 'black', "linewidth", 3, "LineStyle", "--")
% hold on 
% plot(auxdata.time, auxdata.d_min*V_star(:, 1), "color", 'black', "linewidth", 1, "LineStyle", "--")
% hold on 
% plot(auxdata.time, auxdata.d_max*V_star(:, 1) , "color", 'black', "linewidth", 1, "LineStyle", "--")
title("Vehicle -> Leader Headway", "fontsize", 12)
xlabel("Time (s)", "fontsize", 12)
ylabel("Headway (m)", "fontsize", 12)
legend(legend_hw)

%% Calculate Optimized AV Expenditure
% Calculate System Acceleration
optim_system_accel = sum(trapz(auxdata.time, A_star.^2))
% Calculate System Energy Expense
optim_system_energy = sum(trapz(auxdata.time, simplified_fuel_model(V_star, A_star, 'RAV4')))

%%