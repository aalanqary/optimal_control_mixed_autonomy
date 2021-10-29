function l_partial = L_partial(t, i, X, V, U, var, scenario, params)
    global fun_interp
    des_v = params("des_v");
    des_v = des_v(t);
    switch var
        case "x_h"
            if ismember(i+1, scenario("I_a"))
              x = fun_interp(X(:, i), t);
              xf = fun_interp(X(:, i+1), t);
              headway = x - xf - params("l") - params("safe_dist");
              l_partial = (headway < 1) * 2 * log(headway)./headway;
            else
              l_partial = 0;  
            end 
            
        
        % we need x, xl, maybe xf
        case "x_a"
            x = fun_interp(X(:, i), t);
            if i == 1
                xl = scenario("x_leader");
                xl = xl(t);
            else 
                xl = fun_interp(X(:, i-1), t);
            end 
            headway = xl - x - params("l") - params("safe_dist");
            l_partial = (headway < 1) * -2 * log(headway)./headway;
            if headway < 0
                display(i, 'car')
                display(headway, 'headway')
                display(t, 'time')
            end 
            
            if ismember(i+1, scenario("I_a"))
                xf = fun_interp(X(:, i+1), t);
                headway = x - xf - params("l") - params("safe_dist");
                l_partial = l_partial + (headway < 1) * 2 * log(headway)./headway;
            end 
 
        case "v_h"
            v = fun_interp(V(:, scenario("I_h")==i), t);
            l_partial = 2*(v - des_v);
            
        case "v_a"
            l_partial = 0;
%             v = fun_interp(U(:, scenario("I_a")==i), t);
%             l_partial = 2*(v - params("des_v"));
    end 
end 