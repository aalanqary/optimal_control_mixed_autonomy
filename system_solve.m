function [X, V] = system_solve(U, params, scenario)
    global fun_interp forward_sys_solve
%     forward_sys_solve = forward_sys_solve + 1
    U_fun = @(t) fun_interp(U, t);
    [~, XV] = ode15s(@(t,XV) F(t, XV, U_fun(t)', scenario, params), params("t_int"), [scenario("x_0"); scenario("v_0")], params("ode_opt"));
    V = XV(:, end - length(scenario("config")) + 1: end);
    X = XV(:, 1:end - length(scenario("config")));
end 