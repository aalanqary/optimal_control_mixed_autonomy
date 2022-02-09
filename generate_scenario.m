function scenario = generate_scenario(config, x_0, v_0, v_leader, params)
    scenario = containers.Map;
    
    % cars configuration 
    scenario("config") = config;
    scenario("I_a") = find(scenario("config"));
    scenario("I_h") = find(scenario("config") - 1);
    
    % initial conditions 
    scenario("x_0") = x_0(2:end);
    scenario("v_0") = v_0;
    
    % leader velocity 
    scenario("v_leader") = v_leader;
    [~, x_leader] = ode45(@(t, dd) v_leader(t), params("t_int"), x_0(1));
    scenario("x_leader") = @(t) interp1(params("t_int"), x_leader, t);

end 