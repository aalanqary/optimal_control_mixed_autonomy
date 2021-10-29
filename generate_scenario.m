function scenario = generate_scenario(config, x_0, v_0, x_leader, v_leader)
    scenario = containers.Map;
    
    % cars configuration 
    scenario("config") = config;
    scenario("I_a") = find(scenario("config"));
    scenario("I_h") = find(scenario("config") - 1);
    
    % initial conditions 
    if length(x_0) > 1
        scenario("x_0") = x_0(2:end);
    else
        scenario("x_0") = flip(0:length(config) - 1);
        scenario("x_0") = x_0 * scenario("x_0")';
    end 
    
    if length(v_0) > 1
        scenario("v_0") = v_0;
    else
        scenario("v_0") = v_0*ones(length(scenario("I_h")),1);
    end 

    scenario("v_leader") = v_leader;
    scenario("x_leader") = x_leader;

end 