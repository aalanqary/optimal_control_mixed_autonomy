function [X, V] = system_solve(U, params, scenario, is_plot, name)
    global fun_interp
    U_fun = @(t) fun_interp(U, t);
    [~, XV] = ode45(@(t,XV) F(t, XV, U_fun(t), scenario, params), params("t_int"), [scenario("x_0"); scenario("v_0")]);
    V = XV(:, end - length(scenario("I_h")) + 1: end);
    X = XV(:, 1:end - length(scenario("I_h")));
    if is_plot
        plot_x_v(X, V, U, params, scenario, name)
    end 
end 