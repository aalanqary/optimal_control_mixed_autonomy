%% Define problem 
auxdata.platoon = [1, zeros(1, 10)];
auxdata.len_platoon = length(auxdata.platoon);
auxdata.Ia = find(auxdata.platoon);
auxdata.Ih = find(auxdata.platoon - 1);

% Bando-FtL params
auxdata.safe_dist = 2.5; 
auxdata.v_max = 35; 
auxdata.alpha = 0.1;
auxdata.beta = 21;
auxdata.l = 5; 

% objective function auxdata 
auxdata.mu1 = 1;
auxdata.mu2 = 1;

% Constraints auxdata
auxdata.d_min = auxdata.safe_dist;
auxdata.d_max = 50;

%Arctan Barrier auxdata (a(-arctan(bx+c)+pi/2)
auxdata.a = 100;
auxdata.b = 10;
auxdata.c = 1;

T = 250; 
auxdata.utime = (0:1:T)';
auxdata.time = (0:0.1:T)';

% Leader's trajectory short
% vl = @(t) (t<=80) .* 30 ...
%             + (((t>80) + (t<= 120)) ==2) .* (-t./8 + 40) ...
%             + (((t>120) + (t<= 150)) ==2) .* 25 ...
%             + (((t>150) + (t<= 190)) ==2) .* (t./8 + 25/4) ...
%             + (t>190) .* 30;
% vl = vl(auxdata.time); 
% vl = smoothdata(vl, "movmean", 50);
% auxdata.vl = griddedInterpolant(auxdata.time, vl);

% auxdata.vl = @(t) sin(0.1*t) + cos(0.05*t)  - t/100 + 32;



x0 = (eq_headway(auxdata.vl(0), auxdata)+auxdata.l) * flip(0:1:auxdata.len_platoon)';
auxdata.v0 = ones(auxdata.len_platoon, 1) * auxdata.vl(0); 
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[~, xl] = ode45(@(t,x) auxdata.vl(t), auxdata.time, x0(1), opts);
auxdata.x0 = x0(2:end);
auxdata.xl = griddedInterpolant(auxdata.time, xl);

%Initial solution  
U0 = diff(auxdata.vl(auxdata.utime));
U0 = [U0;0];
[X0, V0, A0] = system_solve(U0, auxdata);

%% Run Optimizaer 

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
[U_star, f_val, ~, output, ~, grad] = fmincon(fun, U_star, A, b, Aeq, beq, a_min, a_max, nonlcon, options);
[X_star, V_star, A_star] = system_solve(U_star, auxdata);
toc

