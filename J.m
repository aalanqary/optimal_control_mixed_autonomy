function j = J(X, V, U, scenario, params)
    E = fuel_consumption(X, V, U, scenario, params);
    JT = X(end, :);
    j = params("T")/params("nt") * params("mu1") * sum(E, "all") ...
        - params("mu2") * sum(JT);
end 