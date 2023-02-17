function [time, vl, xl, x0, v0] = generate_scenario(spec)
    time = 0:spec.ts:spec.horizon;
    v1 = 35; v2 = 10; 
    if spec.leader == "wave"
        vl = @(t) sin(0.01 * t) + cos(0.1 *t) + v1;
    elseif spec.leader == "congestion" % 
        t1 = spec.horizon/2; t2 = 3*spec.horizon/4;
        m1 = (v1-v2) / (t1 - t2); 
        vl = @(t) v1 + ...
                  ((t>t1) & (t<= t2)) .* m1 .* (t - t1) + ...
                  (t>t2) .* -(v1-v2);
    elseif spec.leader == "freeflow"
        vl = @(t) v1 + 0.*t; 
    else 
        error("Leader scenario is not defined. Avaliable options: wave, congestion, freeflow")
    end 
    vl(0)
    x0 = spec.th * vl(0) * flip(0:1:spec.len_platoon)'; 
    v0 = ones(spec.len_platoon, 1) * vl(0); 
    opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
    [~, xl] = ode45(@(t,x) vl(t), time, x0(1), opts);
    x0 = x0(2:end);
    xl = griddedInterpolant(time, xl);
end 
