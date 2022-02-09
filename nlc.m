function [c, ceq] = nlc(U, params, scenario)
    [X,V] = system_solve(U, params, scenario);
    xl = scenario("x_leader");
    xl = xl(params("t_int"));
    X_l = [xl, X]; 
    head_way = X_l(:, scenario("I_a")) - X(:, scenario("I_a")) - params("l");
    c_x = sum((min(0, head_way)).^2);
    c_v = sum((min(0, V(:, scenario("I_a")))).^2);
    c = [c_x; c_v];
    ceq = [];
end 