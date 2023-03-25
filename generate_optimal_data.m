auxdata.platoon = [1];
auxdata.len_platoon = length(auxdata.platoon);
auxdata.Ia = find(auxdata.platoon);
auxdata.Ih = find(auxdata.platoon - 1);
auxdata.l = 5;
% Bando-FtL params
auxdata.safe_dist = 4; 
auxdata.v_max = 35; 
auxdata.alpha = 0.1;
auxdata.beta = 21*25;
auxdata.k = 0.2;
auxdata.l = 5; 
% objective function auxdata 
auxdata.mu1 = 1;
auxdata.mu2 = 1;
% Constraints auxdata
auxdata.d_min = 0.5;
auxdata.d_max = 50.0;
files = dir("data");
for i = 3:length(files)
    data = readtable("data/"+files(i).name);
    time = data.Time;
    time = time - time(1);
    T = floor(time(end));
    auxdata.dt = 0.1;
    auxdata.udt = 1; 
    auxdata.utime = (0:auxdata.udt:T)';
    auxdata.time = (0:auxdata.dt:T)';

    vl = data.Velocity * (1000/3600); 
    vl = smoothdata(vl,'movmean',200);
    vl = vl(1:length(auxdata.time));
    auxdata.vl = griddedInterpolant(auxdata.time,vl);
    x0 = (eq_headway(auxdata.vl(0), auxdata)+auxdata.l) * flip(0:1:auxdata.len_platoon)';
    x0 = x0 - x0(2); 
    auxdata.v0 = ones(auxdata.len_platoon, 1) * auxdata.vl(0); 
    opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
    [~, xl] = ode45(@(t,x) auxdata.vl(t), auxdata.time, x0(1), opts);
    auxdata.x0 = x0(2:end);
    auxdata.xl = griddedInterpolant(auxdata.time, xl);
    %Initial solution  
    U0 = diff(auxdata.vl(auxdata.utime)) ./ auxdata.udt;
    U0 = [U0;0];
    [X0, V0, A0] = system_solve(U0, auxdata);

    options = optimoptions('fmincon','Display','iter-detailed', ...
                        'SpecifyObjectiveGradient', true ,...
                        'FunValCheck','on', 'DerivativeCheck', 'off',...
                        'maxfunevals',1e16, 'StepTolerance',1e-12, ...
                        'algorithm', 'sqp');

    fun = @(U) objective_gradient_acc(U, auxdata); 
    nonlcon = []; %@(U) nlc(U, auxdata);
    A = [];
    b = [];
    [A, b] = lc(auxdata); 
    Aeq = []; beq = []; 
    a_min = []; %-3 * ones(size(U0)); 
    a_max = []; %3 * ones(size(U0));
    tic
    [U_star, f_val, ~, output, ~, grad] = fmincon(fun, U0, A, b, Aeq, beq, a_min, a_max, nonlcon, options);
    [X_star, V_star, A_star] = system_solve(U_star, auxdata);
    csvwrite(sprintf('optimal_data/%s', files(i).name),[auxdata.xl(auxdata.time),X_star, vl, V_star, A_star]);
    plot(vl)
    hold on 
    figure(1)
    plot(V_star)
    drawnow;
    toc
end 

