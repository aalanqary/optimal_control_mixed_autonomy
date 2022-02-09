function PQ_dot = F_adjoint(t, PQ, X, V, U, scenario, params)
    global fun_interp
    P = PQ(1:length(scenario("config")));
    Q = PQ(length(scenario("config")) + 1: end);
    P_dot = zeros(length(P), 1);
    Q_dot = zeros(length(Q), 1);
        
    for i = scenario("I_h")
        x = fun_interp(X(:, i), t);
        v = fun_interp(V(:, i), t);
        u = 0; 
        if i == 1 
            xl = scenario("x_leader");
            xl = xl(t);
            vl = scenario("v_leader");
            vl = vl(t); 
        else 
            xl = fun_interp(X(:, i-1), t);
            vl = fun_interp(V(:, i-1), t);
        end
        if ismember(i + 1, scenario("I_h"))
            xf = fun_interp(X(:, i+1), t);
            vf = fun_interp(V(:, i+1), t);
        else
            xf = -1; 
            vf = -1;
        end
        P_dot(i) = L_partial([xl, x, xf], [vl, v, vf], u, "x_h", params) ...
                   - ACC_partial(xl, x, vl, v, 2, params) * Q(i);
        Q_dot(i) = L_partial([xl, x, xf], [vl, v, vf], u, "v_h", params) ...
                                    - P(i) - ACC_partial(xl, x, vl, v, 4, params) * Q(i);
        if ismember(i + 1, scenario("I_h"))
            P_dot(i) = P_dot(i) - ACC_partial(x, xf, v, vf, 1, params) * Q(i+1);
            Q_dot(i) = Q_dot(i) - ACC_partial(x, xf, v, vf, 3, params) * Q(i+1);        
        end 
    end 

    % For x_a need to know: maybe x,v and xf,vf all as functions
    for i = scenario("I_a") 
        x = fun_interp(X(:, i), t);
        v = fun_interp(V(:, i), t);
        u = fun_interp(U(:, scenario("I_a") == i), t);
        xl = -1;
        vl = -1;
        if ismember(i + 1, scenario("I_h"))
            xf = fun_interp(X(:, i+1), t);
            vf = fun_interp(V(:, i+1), t);
        else
            xf = -1; 
            vf = -1;
        end
        P_dot(i) = L_partial([xl, x, xf], [vl, v, vf], u, "x_a", params);
        Q_dot(i) = L_partial([xl, x, xf], [vl, v, vf], u, "v_a", params) - P(i);
        if ismember(i + 1, scenario("I_h"))
            P_dot(i) = P_dot(i) - ACC_partial(x, xf, v, vf, 1, params) * Q(i+1);
            Q_dot(i) = Q_dot(i) - ACC_partial(x, xf, v, vf, 3, params) * Q(i+1);        
        end 
    end
                                
    PQ_dot = [P_dot; Q_dot];
end 