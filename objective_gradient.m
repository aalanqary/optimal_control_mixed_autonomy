function [z, dz] = objective_gradient(U_vec, params, scenario)
    
    % Objective 
    Fu = griddedInterpolant(scenario("time"),U_vec);
    U = @(t) Fu(t);
    [X, V, A] = system_solve(U, params, scenario);
    z = J(X, V, A, scenario, params);
    
    % Gradient
    if nargout >= 2
        Fx = griddedInterpolant(scenario("time"),X);
        Fv = griddedInterpolant(scenario("time"),V);
        [P, Q] = adj_system_solve(@(t) Fx(t), @(t) Fv(t), U, params, scenario);
        dz = -1* Q(:, scenario("Ia")) +  L_partial(X(:, scenario("Ia")), V(:, scenario("Ia")), U(scenario("time")), "u", params);
    end 
 
end 

%%%%% Objective Function %%%%%
function j = J(X, V, A, scenario, params)
    E = simplified_fuel_model(V,A,'RAV4');
    JT = 0; %X(end, :);
    time = scenario("time");
    del_t = time(2) - time(1); 
    j = del_t * params("mu1") * sum(E, "all")/sum(X(end, :), "all") ...
        + params("mu2") * sum(JT);
end 

%%%%% ODE solver for the adjoint system %%%%%
function [P, Q] = adj_system_solve(X, V, U, params, scenario)
    PQ_0 = get_adjoint_ic(X, V, U, params, scenario);
    [~,PQ] = ode15s(@(t,PQ) F_adjoint(t, PQ, X(t), V(t), U(t), scenario, params), flip(scenario("time")), PQ_0, params("ode_opt"));
    PQ = flip(PQ,1);
    P = PQ(:, 1:length(scenario("config")));
    Q = PQ(:, length(scenario("config")) + 1: end);
end 

%%%%% Adjoint system %%%%%
function PQ_dot = F_adjoint(t, PQ, X, V, U, scenario, params)
    P = PQ(1:length(scenario("config")));
    Q = PQ(length(scenario("config")) + 1: end);
    
    P_dot = zeros(length(P), 1);
    Q_dot = zeros(length(Q), 1);
        
    for i = scenario("Ih")
        x = X(i);
        v = V(i);
        u = 0; 
        if i == 1 
            xl = scenario("xl");
            xl = xl(t);
            vl = scenario("vl");
            vl = vl(t); 
        else 
            xl = X(i-1);
            vl = V(i-1);
        end
        if ismember(i + 1, scenario("Ih"))
            xf = X(i+1);
            vf = V(i+1);
        else
            xf = -1; 
            vf = -1;
        end
        P_dot(i) = L_partial([xl, x, xf], [vl, v, vf], u, "x_h", params) ...
                   - ACC_partial(xl, x, vl, v, 2, params) * Q(i);
        Q_dot(i) = L_partial([xl, x, xf], [vl, v, vf], u, "v_h", params) ...
                                    - P(i) - ACC_partial(xl, x, vl, v, 4, params) * Q(i);
        if ismember(i + 1, scenario("Ih"))
            P_dot(i) = P_dot(i) - ACC_partial(x, xf, v, vf, 1, params) * Q(i+1);
            Q_dot(i) = Q_dot(i) - ACC_partial(x, xf, v, vf, 3, params) * Q(i+1);        
        end 
    end 

    % For x_a need to know: maybe x,v and xf,vf all as functions
    for i = scenario("Ia") 
        x = X(i);
        v = V(i);
        u = U(:, scenario("Ia") == i);
        xl = -1;
        vl = -1;
        if ismember(i + 1, scenario("Ih"))
            xf = X(i+1);
            vf = V(i+1);
        else
            xf = -1; 
            vf = -1;
        end
        P_dot(i) = L_partial([xl, x, xf], [vl, v, vf], u, "x_a", params);
        Q_dot(i) = L_partial([xl, x, xf], [vl, v, vf], u, "v_a", params) - P(i);
        if ismember(i + 1, scenario("Ih"))
            P_dot(i) = P_dot(i) - ACC_partial(x, xf, v, vf, 1, params) * Q(i+1);
            Q_dot(i) = Q_dot(i) - ACC_partial(x, xf, v, vf, 3, params) * Q(i+1);        
        end 
    end
    
                                
    PQ_dot = [P_dot; Q_dot];
end 

%%%%% Partial derivatives of the running cost %%%%%
function l_partial = L_partial(X, V, u, var, params)
    load(['RAV4_coeffs.mat']);
    C2 = double(C2); q0 = double(q0);
    
    switch var
        case "x_h" 
            a = ACC(X(1), X(2), V(1), V(2), params);
            partial_a_x = ACC_partial(X(1), X(2), V(1), V(2), 2, params);
            l_partial = (p0 + p1.* V(2) + p2 * V(2).^2) * partial_a_x;
            if a > 0 
                l_partial = l_partial + 2* (q0 + q1 * V(2)) * a * partial_a_x;
            end 
            if V(3) >= 0 %% That is the follower is human car
                af = ACC(X(2), X(3), V(2), V(3), params);
                partial_a_x = ACC_partial(X(2), X(3), V(2), V(3), 1, params);
                l_partial = l_partial + (p0 + p1.* V(3) + p2 * V(3).^2) * partial_a_x;
                if af > 0
                    l_partial = l_partial + 2* (q0 + q1 * V(3)) * af * partial_a_x;
                end 
            end 
            
            
        case "v_h" 
            a = ACC(X(1), X(2), V(1), V(2), params);
            partial_a_v = ACC_partial(X(1), X(2), V(1), V(2), 4, params);
            l_partial = C1 + 2*C2*V(2) + 3*C3*V(2)^2 ...
                        + (p0 + p1.* V(2) + p2 * V(2).^2) * partial_a_v ...
                        + (p1 + 2*p2*V(2)) * a;
            if a > 0 
                l_partial = l_partial + 2* (q0 + q1*V(2)) * a * partial_a_v ...
                            + q1 * a^2;
            end 
            if V(3) >= 0 %% That is the follower is human car
                af = ACC(X(2), X(3), V(2), V(3), params);
                partial_a_v = ACC_partial(X(2), X(3), V(2), V(3), 3, params);
                l_partial = l_partial + (p0 + p1.* V(3) + p2 * V(3).^2) * partial_a_v;
                if af > 0
                    l_partial = l_partial + 2* (q0 + q1 * V(3)) * af * partial_a_v;
                end 
            end 
            
            
        case "x_a" 
            l_partial = 0;
            if V(3) >= 0 %% That is the follower is human car
                af = ACC(X(2), X(3), V(2), V(3), params);
                partial_a_x = ACC_partial(X(2), X(3), V(2), V(3), 1, params);
                l_partial = l_partial + (p0 + p1.* V(3) + p2 * V(3).^2) * partial_a_x;
                if af > 0
                    l_partial = l_partial + 2* (q0 + q1 * V(3)) * af * partial_a_x;
                end 
            end 
            
            
        case "v_a" 
            a = u;
            l_partial = C1 + 2*C2*V(2) + 3*C3*V(2)^2 ...
                        + (p1 + 2*p2*V(2)) * a;
            if a > 0 
                l_partial = l_partial + q1 * a^2;
            end 
            if V(3) >= 0 %% That is the follower is human car
                af = ACC(X(2), X(3), V(2), V(3), params);
                partial_a_v = ACC_partial(X(2), X(3), V(2), V(3), 3, params);
                l_partial = l_partial + (p0 + p1.* V(3) + p2 * V(3).^2) * partial_a_v;
                if af > 0
                    l_partial = l_partial + 2* (q0 + q1 * V(3)) * af * partial_a_v;
                end 
            end 
        
       case "u" 
           l_partial = p0 + p1 .* V + p2 .* V.^2; 
           if u > 0
               l_partial = l_partial + 2.*q0.*u + 2.*q1.*V.*u;
           end 
                     
    end 
end 

%%%%% Partial derivatives of the ACC function %%%%%
function a_partial = ACC_partial(xl, xf, vl, vf, var, params)
    
    V_prime = @(head_way) (params("v_max") * (sech(head_way - params("safe_dist"))).^2) ./ (1+tanh(params("l")+params("safe_dist")));
    head_way = xl - xf - params("l");
    switch var
        case 1 
            a_partial = params("alpha") * V_prime(head_way) - 2*params("beta")*((vl-vf)./head_way.^3);
        case 2
            a_partial = -1*params("alpha") * V_prime(head_way) + 2*params("beta")*((vl-vf)./head_way.^3);
        case 3
            a_partial = params("beta") ./ head_way.^2;
        case 4
            a_partial = -1*params("alpha") - params("beta") ./ head_way.^2;
    end 
end 

function PQ_0 = get_adjoint_ic(X, V, U, params, scenario)
    Q_0 = zeros(length(scenario("config")), 1);
    P_0 = zeros(length(scenario("config")), 1);
    PQ_0 = [P_0;Q_0];
end 