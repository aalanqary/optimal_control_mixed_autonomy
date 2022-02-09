function PQ_0 = get_adjoint_ic(X, V, U, params, scenario)
    Q_0 = zeros(length(scenario("config")), 1);
    P_0 = ones(length(scenario("config")), 1);
    PQ_0 = [P_0;Q_0];
end 