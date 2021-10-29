function plot_x_v(X, V, U, params, scenario, name)
    f = figure('Color','none', 'Unit', 'inches', 'Position', [0,0,15,5]);
    labels = ["Leader"];
    hcar = 1;
    acar = 1;
    for car = scenario("config")
        if car
            labels = [labels, "Autonomous " + acar];
            acar = acar + 1;
        else
            labels = [labels, "Human " + hcar];
            hcar = hcar + 1;
        end 
    end
    
    subplot(1,2,1)
    x_leader = scenario("x_leader");
    plot(params("t_int"),x_leader(params("t_int")),params("t_int"),X)
    legend(labels)
    xlabel("time (s)", 'FontSize', 15)
    xlabel("Position (m)", 'FontSize', 15)
    title("Positions", 'FontSize', 20)
    
    subplot(1,2,2)
    v_leader = scenario("v_leader");
    V_ordered = zeros(params("nt"), length(scenario("config")));
    V_ordered(:, scenario("I_a")) = U;
    V_ordered(:, scenario("I_h")) = V;
    plot(params("t_int"),v_leader(params("t_int")),params("t_int"),V_ordered)
    legend(labels)
    xlabel("time (s)", 'FontSize', 15)
    xlabel("Velocity (m/s)", 'FontSize', 15)
    title("Velocities", 'FontSize', 20)
    
    saveas(f, sprintf('scenarios/%s_%i.png', name, scenario("number")))

end 