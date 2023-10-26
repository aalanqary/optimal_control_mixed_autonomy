function [auxdata, leader] = problem_auxdata(platoon, constraints, leader_type)
%     platoon: list of 1 and 0 representing locations of the AV and HV
%     constraints: type of problem in terms of dealing with constraints
%     leader_type: the name of the leader trajectory ["simple", "sin", "real"]

    %Platoon configeration
    auxdata.platoon = platoon;
    auxdata.len_platoon = length(platoon);
    auxdata.Ia = find(platoon);
    auxdata.Ih = find(platoon - 1);
    
    %Bando-FtL params
    auxdata.safe_dist = 4; 
    auxdata.v_max = 35; 
    auxdata.alpha = 0.1;
    auxdata.beta = 21*25;
    auxdata.k = 0.2;
    auxdata.l = 5; 
    
    %objective function 

    %constraints 
    if constraints == "linear"
        auxdata.h_min = 0.5; %time headway  
        auxdata.h_max = 3; %time headway 
    elseif constraints == "penalty_minmax" | constraints == "greedy"
        auxdata.d_min = 5;
        auxdata.d_max = 120;
        auxdata.mu_min = 1e-1; 
        auxdata.mu_max = 1e-1; 
    elseif constraints == "smooth_penalty"
        auxdata.d_min = 5;
        auxdata.eps = 0.0; 
        auxdata.mu_min = 0.0001; 
    else 
        error("wrong consstraints choice")
    end 
    %Leader trajectory
    if leader_type == "simple" 
        T = 900; 
        auxdata.dt = 0.1;
        auxdata.udt = 1; 
        auxdata.utime = (0:auxdata.udt:T)';
        auxdata.time = (0:auxdata.dt:T)';
        
        vl = @(t) (t<=120) .* 30 ...
                    + (((t>120) + (t<= 300)) ==2) .* (-t./9 + 130/3) ...
                    + (((t>300) + (t<= 600)) ==2) .* 10 ...
                    + (((t>600) + (t<= 780)) ==2) .* (t./9 - 170/3) ...
                    + (t>780) .* 30;
        vl = vl(auxdata.time); 
        vl = smoothdata(vl, "movmean", 400);
        leader.v = griddedInterpolant(auxdata.time, vl);
    elseif leader_type == "simplest"
        T = 10; 
        auxdata.dt = 0.1;
        auxdata.udt = 1; 
        auxdata.utime = (0:auxdata.udt:T)';
        auxdata.time = (0:auxdata.dt:T)';
        leader.v = @(t) -0.5*t + 30;
    elseif leader_type == "sin"
        T = 250; 
        auxdata.dt = 0.1;
        auxdata.udt = 1; 
        auxdata.utime = (0:auxdata.udt:T)';
        auxdata.time = (0:auxdata.dt:T)';
        
        leader.v = @(t) sin(0.1*t) + cos(0.05*t)  - t/100 + 32;
    else 
        T = 704; 
        auxdata.dt = 0.5;
        auxdata.udt = 5; 
        auxdata.utime = (0:auxdata.udt:T)';
        auxdata.time = (0:auxdata.dt:T)';
        
        data = readtable("data_v2_preprocessed_west/2021-04-22-12-47-13_2T3MWRFVXLW056972_masterArray_0_7050.csv");
        vl = data.Velocity * (1000/3600); 
        vl = smoothdata(vl,'movmean',200);
        leader.v = griddedInterpolant(data.Time - min(data.Time),vl);
    end 
    eq = round(eq_headway(leader.v(0), auxdata), 5);
    x0 = (eq+auxdata.l) * flip(0:1:auxdata.len_platoon)';
    auxdata.v0 = ones(auxdata.len_platoon, 1) * leader.v(0); 
    opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
    [~, xl] = ode45(@(t,x) leader.v(t), auxdata.time, 0, opts);
    xl = xl + x0(1);
    auxdata.x0 = x0(2:end);
    leader.x = griddedInterpolant(auxdata.time, xl);
end 

function h = eq_headway(v, auxdata)
    C = tanh(auxdata.l + auxdata.safe_dist);
    h = atanh((v * (1+C)./auxdata.v_max) - C) + auxdata.safe_dist;
    h = h/auxdata.k;
end 