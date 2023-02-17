function scenario = prepare_data(x0_rule, platoon, nt)
    scenario = containers.Map;
    scenario("config") = platoon;
    scenario("Ia") = find(scenario("config"));
    scenario("Ih") = find(scenario("config") - 1);
    len_platoon = length(platoon);
    
    % Synthetic data 
    if true 
        time = 0:0.1:500;
        vl_fun = @(t) (t < 300) .* 35 + (t >= 300) .* (35-(t-300)/10);
    end 
    % Read Data
    if false
        data = readtable(sprintf("%s/%s", data_dir, file));
        
        time = data.Time;
        if nargin < 5
            nt = length(time);
        end 
        time = time(1:nt) - time(1);
        
        vl = data.Velocity; 
        vl = vl(1:nt) * (1000/3600);
        Fv = griddedInterpolant(time,vl);
        vl_fun = @(t) Fv(t);
    end 

    x0 = x0_rule * vl_fun(time(1)) * flip(0:1:len_platoon)'; 
    v0 = ones(len_platoon, 1) * vl_fun(time(1)); 
    
    opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
    [~, xl] = ode45(@(t,x) vl_fun(t), time, x0(1), opts);
    Fx = griddedInterpolant(time, xl);
    xl_fun = @(t) Fx(t);  
    
    scenario("time") = time;
    scenario("x0") = x0(2:end);
    scenario("v0") = v0;
    scenario("vl") = vl_fun;
    scenario("xl") = xl_fun;
end 