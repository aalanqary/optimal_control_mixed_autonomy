function [P, Q] = back_system_solve(X, V, U, params, scenario)
    PQ_0 = get_adjoint_ic(X, V, U, params, scenario);
    [~,PQ] = ode15s(@(t,PQ) F_adjoint(t, PQ, X, V, U, scenario, params), flip(params("t_int")), PQ_0, params("ode_opt"));
    PQ = flip(PQ,1);
    P = PQ(:, 1:length(scenario("config")));
    Q = PQ(:, length(scenario("config")) + 1: end);
end 