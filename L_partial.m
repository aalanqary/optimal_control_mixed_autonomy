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