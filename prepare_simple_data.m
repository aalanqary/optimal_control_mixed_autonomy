function scenario_simple = prepare_simple_data(data_dir, file, x0_rule, platoon, nt)
    % change so that the inital x0 = v0=v(end)=0 
    scenario_simple = containers.Map;
    scenario_simple("config") = platoon;
    scenario_simple("Ia") = find(scenario_simple("config"));
    scenario_simple("Ih") = find(scenario_simple("config") - 1);
    len_platoon = length(platoon);
    
    % Read Data
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
   
    x0 = x0_rule * vl(1) * flip(0:1:len_platoon)'; 
    v0 = ones(len_platoon, 1) * vl(1); 
    
    opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
    [~, xl] = ode45(@(t,x) vl_fun(t), time, x0(1), opts);
    Fx = griddedInterpolant(time, xl);
    xl_fun = @(t) Fx(t);  
    
    scenario_simple("time") = time;
    scenario_simple("x0") = x0(2:end);
    
    scenario_simple("v0") = v0;
    scenario_simple("vl") = vl_fun;
    scenario_simple("xl") = xl_fun;
end 