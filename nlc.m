function [c, ceq] = nlc(U_vec, params, scenario)
    Fu = griddedInterpolant(scenario("time"),U_vec);
    U = @(t) Fu(t);
    [X,V] = system_solve(U, params, scenario);
    xl = scenario("xl");
    xl = xl(scenario("time"));
    X_l = [xl, X]; 
    head_way = X_l(:, scenario("Ia")) - X(:, scenario("Ia")) - params("l");
    c_x = sum((min(0, head_way)).^2);
    c_v = sum((min(0, V(:, scenario("Ia")))).^2);
    c = [c_x; c_v];
    ceq = [];
end 