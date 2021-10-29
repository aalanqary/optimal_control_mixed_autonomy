function PQ_dot = F_adjoint(t, PQ, X, V, U, scenario, params) %X\in R(t_int X N+M), V \in R(t_int X M), U \in R(t_int X N)
    global fun_interp
    Q = PQ(end - length(scenario("I_h")) + 1: end);
    P = PQ(1:end - length(scenario("I_h")));
    
    P_dot = zeros(length(P), 1);
    Q_dot = zeros(length(Q), 1);
        
    % For x_h and v_h we need to know: x,v, xl, vl maybe xf,vf all as functions
    for i = scenario("I_h")
        x = fun_interp(X(:, i), t);
        v = fun_interp(V(:, scenario("I_h")==i), t);
        if i == 1 
            xl = scenario("x_leader");
            xl = xl(t);
            vl = scenario("v_leader");
            vl = vl(t); 
        elseif ismember(i - 1, scenario("I_h"))
            xl = fun_interp(X(:, i-1), t);
            vl = fun_interp(V(:, scenario("I_h")==i-1), t);
        else 
            xl = fun_interp(X(:, i-1), t);
            vl = fun_interp(U(:, scenario("I_a")==i-1), t);
        end
        P_dot(i) = L_partial(t, i, X, V, U, "x_h", scenario, params) ...
                   - ACC_partial(xl, x, vl, v, 2, params) * Q(scenario("I_h")==i);
        Q_dot(scenario("I_h")==i) = L_partial(t, i, X, V, U, "v_h", scenario, params) ...
                                    - P(i) ...
                                    - ACC_partial(xl, x, vl, v, 4, params) * Q(scenario("I_h")==i);
        if ismember(i + 1, scenario("I_h"))
            xf = fun_interp(X(:, i+1), t);
            vf = fun_interp(V(:, scenario("I_h")==i+1), t);
            P_dot(i) = P_dot(i) - ACC_partial(x, xf, v, vf, 1, params) * Q(scenario("I_h")==i+1);
            Q_dot(scenario("I_h")==i) = Q_dot(scenario("I_h")==i) - ACC_partial(x, xf, v, vf, 3, params) * Q(scenario("I_h")==i+1);        
        end 
    end 

    % For x_a need to know: maybe x,v and xf,vf all as functions
    for i = scenario("I_a") 
        P_dot(i) = L_partial(t, i, X, V, U, "x_a", scenario, params);
        if ismember(i+1, scenario("I_h"))
            x = fun_interp(X(:, i), t);
            v = fun_interp(U(:, scenario("I_a")==i), t);
            xf = fun_interp(X(:, i+1), t);
            vf = fun_interp(V(:, scenario("I_h")==i+1), t);
            P_dot(i) = P_dot(i) - ACC_partial(x, xf, v, vf, 1, params) * Q(scenario("I_h")==i+1);
        end 
    end
    PQ_dot = [P_dot; Q_dot];
end 