function [X_0, V_0, U_0] = initial_solution(scenario, params)
    original_I_a = scenario("I_a");
    scenario("config") = zeros(1, length(scenario("config")));
    scenario("I_a") = find(scenario("config"));
    scenario("I_h") = find(scenario("config") - 1);
    params("l") = 2 * params("l");
    [X_0, V_0] = system_solve(zeros(params("nt"), 1), params, scenario);
    xl = scenario("x_leader");
    xl = xl(params("t_int"));
    vl = scenario("v_leader");
    vl = vl(params("t_int"));
    X_l = [xl, X_0]; 
    V_l = [vl, V_0]; 
    U_0 = ACC(X_l(:, original_I_a), X_0(:, original_I_a), ... 
              V_l(:, original_I_a), V_0(:, original_I_a), params);
end 