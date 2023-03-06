%% Define problem 
auxdata.platoon = [1,0,0];
auxdata.len_platoon = length(auxdata.platoon);
auxdata.Ia = find(auxdata.platoon);
auxdata.Ih = find(auxdata.platoon - 1);

% Bando-FtL params
auxdata.safe_dist = 2.5; 
auxdata.v_max = 35; 
auxdata.alpha = 0.5;
auxdata.beta = 21;
auxdata.l = 5; 

% objective function auxdata 
auxdata.mu1 = 2;
auxdata.mu2 = 0.2;
auxdata.iter = 0;

% Constraints auxdata
auxdata.d_min = auxdata.safe_dist;

% optimizer auxdata 
auxdata.eps = 2;
auxdata.gamma = 120;

%Arctan Barrier auxdata (a(-arctan(bx+c)+pi/2)
auxdata.a = 100;
auxdata.b = 10000;
auxdata.c = 100;

% Leader's trajectory long
    % auxdata.vl = @(t) (t<=120) * 30 ...
    %             + (((t>120) + (t<= 240)) ==2) .* (-t/6 + 50) ...
    %             + (((t>240) + (t<= 420)) ==2) .* 10 ...
    %             + (((t>420) + (t<= 540)) ==2) .* (t/6 - 60) ...
    %             + (t>540) .* 30;

% Leader's trajectory short
auxdata.vl = @(t) (t<=80) .* 30 ...
            + (((t>80) + (t<= 120)) ==2) .* (-t./8 + 40) ...
            + (((t>120) + (t<= 150)) ==2) .* 25 ...
            + (((t>150) + (t<= 190)) ==2) .* (t./8 + 25/4) ...
            + (t>190) .* 30;

auxdata.utime = (0:1:240)';
auxdata.time = (0:0.1:240)';
x0 = (eq_headway(auxdata.vl(0), auxdata)+auxdata.l) * flip(0:1:auxdata.len_platoon)';
auxdata.v0 = ones(auxdata.len_platoon, 1) * auxdata.vl(0); 
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[~, xl] = ode45(@(t,x) auxdata.vl(t), auxdata.time, x0(1), opts);
auxdata.x0 = x0(2:end);
auxdata.xl = griddedInterpolant(auxdata.time, xl);

%Initial solution  
U0 = diff(auxdata.vl(auxdata.utime) - 5);
U0 = [U0;0];
[X0, V0, A0] = system_solve(U0, auxdata);

%% Run Optimizer (first iteration)

% First iteration

options = optimoptions('fmincon','Display','iter-detailed', ...
                        'SpecifyObjectiveGradient', true ,...
                        'FunValCheck','on', 'DerivativeCheck', 'off',...
                        'maxfunevals',1e6, 'StepTolerance',1e-12, ...
                        'algorithm', 'sqp', ...
                        'MaxIterations', 400);

fun = @(U) objective_gradient_acc(U, auxdata); 
nonlcon = []; %@(U) nlc(U, auxdata);
A = [];
b = [];
[A, b] = lc(auxdata); 
Aeq = []; beq = []; 
a_min = -3 * ones(size(U0)); 
a_max = 3 * ones(size(U0));
[U_star, f_val, ~, output, ~, grad] = fmincon(fun, U0, A, b, Aeq, beq, a_min, a_max, nonlcon, options);
[X_star, V_star] = system_solve(U_star, auxdata);
figure(5)
plot(U_star)
drawnow;
title("U_star")
figure(6)
Xl = [auxdata.xl(auxdata.time), X_star];
plot(Xl(:, 1) - Xl(:, 2) - auxdata.l)
title("AV Headway")
drawnow;

% %% Second iteration
% options = optimoptions('fmincon','Display','iter-detailed', ...
%                         'SpecifyObjectiveGradient', true ,...
%                         'FunValCheck','on', 'DerivativeCheck', 'off',...
%                         'maxfunevals',1e6, 'StepTolerance',1e-12, ...
%                         'algorithm', 'sqp', ...
%                         'MaxIterations', 500);
% auxdata.mu2 = 1;
% fun = @(U) objective_gradient_acc(U, auxdata); 
% nonlcon = []; %@(U) nlc(U, auxdata);
% A = [];
% b = [];
% [A, b] = lc(auxdata); 
% Aeq = []; beq = []; 
% a_min = -3 * ones(size(U0)); 
% a_max = 3 * ones(size(U0));
% [U_star_fin, f_val, ~, output, ~, grad] = fmincon(fun, U_star, A, b, Aeq, beq, a_min, a_max, nonlcon, options);
% [X_star, V_star] =system_solve(U_star_fin, auxdata); 
% figure(7)
% plot(U_star_fin)
% drawnow;
% title("U_star Second iteration")
% figure(8)
% Xl = [auxdata.xl(auxdata.time), X_star];
% plot(Xl(:, 1) - Xl(:, 2) - auxdata.l)
% title("AV Headway")
% drawnow;
% figure(9)
% Vl = [auxdata.vl(auxdata.time), V]; 
% plot(Vl)
% drawnow;