%% Define problem 

auxdata.platoon = [1, zeros(1, 25)];
auxdata.len_platoon = length(auxdata.platoon);
auxdata.Ia = find(auxdata.platoon);
auxdata.Ih = find(auxdata.platoon - 1);


% Leader's trajectory
% auxdata.utime = (0:1:50)';
% auxdata.time = (0:0.1:50)';
% auxdata.vl = @(t) 30 + 0.*t; 
% x0 = 0.27 * auxdata.vl(0) * flip(0:1:auxdata.len_platoon)'; 
% auxdata.v0 = ones(auxdata.len_platoon, 1) * auxdata.vl(0); 
% opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
% [~, xl] = ode45(@(t,x) auxdata.vl(t), auxdata.time, x0(1), opts);
% auxdata.x0 = x0(2:end);
% auxdata.xl = griddedInterpolant(auxdata.time, xl);

% Bando-FtL params
auxdata.safe_dist = 3; 
auxdata.v_max = 35; 
auxdata.alpha = 0.11;
auxdata.beta = 30;
auxdata.l = 4.5; 

% objective function auxdata 
auxdata.mu1 = 1;
auxdata.mu2 = 1;
auxdata.v_des = 28; 

% Constraints auxdata
auxdata.d_min = auxdata.safe_dist;

% optimizer auxdata 
auxdata.eps = 3;
auxdata.gamma = 120;

% Leader's trajectory
T = 704.9;
auxdata.utime = (0:1:T)';
auxdata.time = (0:0.1:T)';
% auxdata.vl = @(t) (t<=120) * 30 ...
%             + (((t>120) + (t<= 240)) ==2) .* (-t/6 + 50) ...
%             + (((t>240) + (t<= 420)) ==2) .* 10 ...
%             + (((t>420) + (t<= 540)) ==2) .* (t/6 - 60) ...
%             + (t>540) .* 30;

% auxdata.vl = @(t) (t<=100) .* (-t./60 + 30) ...
%             + (((t>120) + (t<= 360)) ==2) .* (-t./30 + 32) ...
%             + (((t>360) + (t<= 480)) ==2) .* 20 ...
%             + (((t>480) + (t<= 840)) ==2) .* (t./36 + 20/3) ...
%             + (t>840) .* 30;

% auxdata.vl = @(t) (t<=80) .* 30 ...
%             + (((t>80) + (t<= 120)) ==2) .* (-t./8 + 40) ...
%             + (((t>120) + (t<= 150)) ==2) .* 25 ...
%             + (((t>150) + (t<= 190)) ==2) .* (t./8 + 25/4) ...
%             + (t>190) .* 30;

data = readtable("2021-04-22-12-47-13_2T3MWRFVXLW056972_masterArray_0_7050.csv");
vl = data.Velocity * (1000/3600); 
% vl = auxdata.vl(auxdata.time);
% vl = smoothdata(vl,'movmean',200);
auxdata.vl = griddedInterpolant(auxdata.time,vl);

x0 = (eq_headway(auxdata.vl(0), auxdata)+auxdata.l) * flip(0:1:auxdata.len_platoon)';
auxdata.v0 = ones(auxdata.len_platoon, 1) * auxdata.vl(0); 
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[~, xl] = ode45(@(t,x) auxdata.vl(t), auxdata.time, x0(1), opts);
auxdata.x0 = x0(2:end);
auxdata.xl = griddedInterpolant(auxdata.time, xl);


%Initial solution  
% U0 = diff(auxdata.vl(auxdata.utime) - 5);
U0 = diff(smoothdata(auxdata.vl(auxdata.utime), 'movmean', 10)) ./ diff(auxdata.utime);
% U0 = diff(auxdata.vl(auxdata.utime)) ./ diff(auxdata.utime);
U0 = [U0;0];
[X0, V0, A0] = system_solve(U0, auxdata);
[X_star, V_star, A_star] = system_solve(U_star, auxdata);
% [X0, V0, U0] = initial_solution(auxdata);
% U0 = griddedInterpolant(auxdata.time, U0);
% U0 = U0(auxdata.utime);
% [X0, V0] = system_solve(U0, auxdata);
% U0 = zeros(length(auxdata.utime), length(auxdata.Ia));
% [X0, V0] = system_solve(U0, auxdata);

%% Run Optimizaer 

options = optimoptions('fmincon','Display','iter-detailed', ...
                        'SpecifyObjectiveGradient', true ,...
                        'FunValCheck','on', 'DerivativeCheck', 'off',...
                        'maxfunevals',1e6, 'StepTolerance',1e-12, ...
                        'algorithm', 'sqp');

fun = @(U) objective_gradient_acc(U, auxdata); 
nonlcon = []; %@(U) nlc(U, auxdata);
A = [];
b = [];
[A, b] = lc(auxdata); 
Aeq = []; beq = []; 
a_min = [];%-3 * ones(size(U0)); 
a_max = []; %3 * ones(size(U0));
[U_star, f_val, ~, output, ~, grad] = fmincon(fun, U0, A, b, Aeq, beq, a_min, a_max, nonlcon, options);
[X_star, V_star, A_star] = system_solve(U_star, auxdata);

