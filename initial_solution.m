function [X, V, U_0] = initial_solution(scenario, params)
options_ode = [];
    original_I_a = scenario("I_a");
    scenario("config") = zeros(1, length(scenario("config")));
    scenario("I_a") = find(scenario("config"));
    scenario("I_h") = find(scenario("config") - 1);
    scenario("v_0") = 21 * ones(length(scenario("config")), 1);
    params("l") = params("l") + 1.1 * params("safe_dist");
    [~, XV] = ode45(@(t,XV) F(t, XV, 0, scenario, params), params("t_int"), [scenario("x_0"); scenario("v_0")], options_ode);
    V = XV(:, end - length(scenario("I_h")) + 1: end);
    X = XV(:, 1:end - length(scenario("I_h")));
    U_0 = V(:, original_I_a);
end 