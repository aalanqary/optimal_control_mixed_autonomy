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




%%Plot initial AV velocity  
vl = auxdata.vl(auxdata.time);
figure()
plot(auxdata.time, vl, "color", '#5E5E5E', "linewidth", 1)
hold on 
plot(auxdata.time, V0(:, 1), "color", '#EE220C', "linewidth", 2, "LineStyle", "--")
title("AV Initialization", "fontsize", 12)
xlabel("Time (s)")
ylabel("Velocity (m/s^2)")

%%Plot initial AV acceleration  
vl = auxdata.vl(auxdata.time);
al = diff(vl)./diff(auxdata.time);
al = [al; al(end)];
figure()
plot(auxdata.time, al, "color", '#5E5E5E', "linewidth", 1)
hold on 
plot(auxdata.time, A0(:, 1), "color", '#EE220C', "linewidth", 2, "LineStyle", "--")
title("AV Initialization", "fontsize", 12)
xlabel("Time (s)")
ylabel("Velocity (m/s^2)")




%%Plot HV simulation (v)
figure()
vl = auxdata.vl(auxdata.time);
al = diff(vl)./diff(auxdata.time);
al = [al; al(end)];
figure()
plot(auxdata.time, V0, "color", '#3274B5', "linewidth", 0.5)
hold on
plot(auxdata.time, vl, "color", '#5E5E5E', "linewidth", 1.5)
title("Simulation of human car-following behaviour", "fontsize", 12)
xlabel("Time (s)", "fontsize", 12)
ylabel("Velocity (m/s)", "fontsize", 12)


%%Plot HV simulation (a)
figure()
vl = auxdata.vl(auxdata.time);
al = diff(vl)./diff(auxdata.time);
al = [al; al(end)];
plot(auxdata.time, A0, "color", '#3274B5', "linewidth", 0.5)
hold on 
plot(auxdata.time, al, "color", '#5E5E5E', "linewidth", 1.5)
title("Simulation of human car-following behaviour", "fontsize", 12)
xlabel("Time (s)", "fontsize", 12)
ylabel("Acceleration (m/s^2)", "fontsize", 12)





%%Plot optimal simulation (v)
figure()
plot(auxdata.time, V_star(:, 2:end), "color", '#3274B5', "linewidth", 0.5)
hold on 
plot(auxdata.time, V_star(:, 1), "color", '#EE220C', "linewidth", 1)
hold on
plot(auxdata.time, auxdata.vl(auxdata.time), "color", '#5E5E5E', "linewidth", 1)
title("Simulation of mixed-autonomy car-following behaviour", "fontsize", 12)
xlabel("Time (s)", "fontsize", 12)
ylabel("Velocity (m/s)", "fontsize", 12)

%%Plot optimal simulation (a)
figure()
plot(auxdata.time, A_star(:, 2:end), "color", '#3274B5', "linewidth", 0.5)
hold on 
plot(auxdata.time, A_star(:, 1), "color", '#EE220C', "linewidth", 1)
hold on
vl = auxdata.vl(auxdata.time);
al = diff(vl)./diff(auxdata.time);
al = [al; al(end)];
plot(auxdata.time, al, "color", '#5E5E5E', "linewidth", 1)
title("Simulation of mixed-autonomy car-following behaviour", "fontsize", 12)
xlabel("Time (s)", "fontsize", 12)
ylabel("Acceleration (m/s)", "fontsize", 12)

%%Plot optimal simulation (headway)
figure()
xl = auxdata.xl(auxdata.time);
plot(auxdata.time, xl - X_star(:, 1) - auxdata.l, "color", '#EE220C', "linewidth", 1)
title("The headway between the AV and the leader", "fontsize", 12)
xlabel("Time (s)", "fontsize", 12)
ylabel("Acceleration (m/s)", "fontsize", 12)




% Plot comparision between optimal results 
figure()
plot(auxdata.utime, UAV + 0.1, "linewidth", 1)
hold on 
plot(auxdata.utime, U10 + 0.05, "linewidth", 1)
hold on 
plot(auxdata.utime, U30, "linewidth", 1)

title("comparision of optimal solutions", "fontsize", 12)
xlabel("Time (s)")
ylabel("Acceleration (m/s^2)")
legend(["0AV", "10AV", "30AV"])
