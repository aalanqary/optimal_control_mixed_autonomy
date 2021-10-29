function l_partial = L_partial(t, i, X, V, U, var, scenario, params)
    global fun_interp
    load(['RAV4_coeffs.mat']);
    des_v = params("des_v");
    des_v = des_v(t);
    switch var
        case "x_h" 
            % Get car states
            x = fun_interp(X(:, i), t);
            v = fun_interp(V(:, scenario("I_h")==i), t);
            % Get leader states
            if i == 1
                xl = scenario("x_leader");
                xl = xl(t);
                vl = scenario("v_leader");
                vl = vl(t); 
            elseif ismember(i-1, scenario("I_a"))
                xl = fun_interp(X(:, i-1), t);
                vl = fun_interp(U(:, scenario("I_a")==i-1), t);
            else 
                xl = fun_interp(X(:, i-1), t);
                vl = fun_interp(V(:, scenario("I_h")==i-1), t);
            end 
            % compute base l_partial 
            a = ACC(xl, x, vl, v, params); 
            l_partial = (p0 + p1 .* v + p2.*v.^2) .* ACC_partial(xl, x, vl, v, 2, params); 
            if a > 0
                l_partial = l_partial + 2.*(q0 + q1 .* v) * a * ACC_partial(xl, x, vl, v, 2, params);
            end 
            % compute additional terms if follower is human or autonomouse 
            if ismember(i+1, scenario("I_h"))
                xf = fun_interp(X(:, i+1), t);
                vf = fun_interp(V(:, scenario("I_h")==i+1), t);
                af = ACC(x, xf, v, vf, params); 
                l_partial = l_partial + (p0 + p1 .* vf + p2.*vf.^2) .* ACC_partial(x, xf, v, vf, 1, params); 
                if af > 0
                    l_partial = l_partial + 2.*(q0 + q1 .* vf) * af * ACC_partial(x, xf, v, vf, 1, params);
                end 
            elseif ismember(i+1, scenario("I_a"))
                xf = fun_interp(X(:, i+1), t);
                headway = x - xf - params("l") - params("safe_dist");
                if headway < 0.1
                    t
                    headway
                end 
                l_partial = l_partial + (headway < 1) * 2 * log(headway)./headway;
            end 
           
            
        case "v_h"
            % Get car states
            x = fun_interp(X(:, i), t);
            v = fun_interp(V(:, scenario("I_h")==i), t);
            % Get leader states
            if i == 1
                xl = scenario("x_leader");
                xl = xl(t);
                vl = scenario("v_leader");
                vl = vl(t); 
            elseif ismember(i-1, scenario("I_a"))
                xl = fun_interp(X(:, i-1), t);
                vl = fun_interp(U(:, scenario("I_a")==i-1), t);
            else 
                xl = fun_interp(X(:, i-1), t);
                vl = fun_interp(V(:, scenario("I_h")==i-1), t);
            end 
            % compute base l_partial 
            a = ACC(xl, x, vl, v, params); 
            l_partial = 2*(v - des_v) ...
                      + C1 + 2 .*C2 .* v ...
                      + (p0 + p1 .* v + p2.*v.^2) .* ACC_partial(xl, x, vl, v, 4, params) + (p1 + 2 .* p2 .*v) * a; 
            if a > 0
                l_partial = l_partial + 2.*(q0 + q1 .* v) * a * ACC_partial(xl, x, vl, v, 4, params) + q1 .* a.^2;
            end 
            % compute additional terms if follower is human or autonomouse 
            if ismember(i+1, scenario("I_h"))
                xf = fun_interp(X(:, i+1), t);
                vf = fun_interp(V(:, scenario("I_h")==i+1), t);
                af = ACC(x, xf, v, vf, params); 
                l_partial = l_partial + (p0 + p1 .* vf + p2.*vf.^2) .* ACC_partial(x, xf, v, vf, 3, params); 
                if af > 0
                    l_partial = l_partial + 2.*(q0 + q1 .* vf) * af * ACC_partial(x, xf, v, vf, 3, params);
                end 
            end 
                       
            
        case "x_a"
            % Get car states
            x = fun_interp(X(:, i), t);
            v = fun_interp(U(:, scenario("I_a")==i), t);
            % Get leader states
            if i == 1
                xl = scenario("x_leader");
                xl = xl(t);
                vl = scenario("v_leader");
                vl = vl(t); 
            elseif ismember(i-1, scenario("I_a"))
                xl = fun_interp(X(:, i-1), t);
                vl = fun_interp(U(:, scenario("I_a")==i-1), t);
            else 
                xl = fun_interp(X(:, i-1), t);
                vl = fun_interp(V(:, scenario("I_h")==i-1), t);
            end 
            % compute base l_partial 
            headway = xl - x - params("l") - params("safe_dist");
            l_partial = (headway < 1) * -2 * log(headway)./headway;
            % compute additional terms if follower is human or autonomouse 
            if ismember(i+1, scenario("I_h"))
                xf = fun_interp(X(:, i+1), t);
                vf = fun_interp(V(:, scenario("I_h")==i+1), t);
                af = ACC(x, xf, v, vf, params); 
                l_partial = l_partial + (p0 + p1 .* vf + p2.*vf.^2) .* ACC_partial(x, xf, v, vf, 1, params); 
                if af > 0
                    l_partial = l_partial + 2.*(q0 + q1 .* vf) * af * ACC_partial(x, xf, v, vf, 1, params);
                end 
            elseif ismember(i+1, scenario("I_a"))
                xf = fun_interp(X(:, i+1), t);
                headway = x - xf - params("l") - params("safe_dist");
                l_partial = l_partial + (headway < 1) * 2 * log(headway)./headway;
            end 
            
            
        case "v_a"
            % compute base l_partial 
            l_partial = 0;
            if ismember(i+1, scenario("I_h"))
                x = fun_interp(X(:, i), t);
                v = fun_interp(U(:, scenario("I_a")==i), t);
                xf = fun_interp(X(:, i+1), t);
                vf = fun_interp(V(:, scenario("I_h")==i+1), t);
                af = ACC(x, xf, v, vf, params); 
                l_partial = l_partial + (p0 + p1 .* vf + p2.*vf.^2) .* ACC_partial(x, xf, v, vf, 3, params); 
                if af > 0
                    l_partial = l_partial + 2.*(q0 + q1 .* vf) * af * ACC_partial(x, xf, v, vf, 3, params);
                end 
            end 
    end 
end 