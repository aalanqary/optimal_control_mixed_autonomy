%% Define problem 
auxdata.platoon = [1,1,1,0];
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
auxdata.mu1 = 2;
auxdata.mu2 = 0.4;
auxdata.iter = 0;

% Constraints auxdata
auxdata.d_min = auxdata.safe_dist;

% optimizer auxdata 
auxdata.eps = 2;
auxdata.gamma = 120;

%Arctan Barrier auxdata (a(-arctan(bx+c)+pi/2)
auxdata.a = 1;
auxdata.b = 10000000;
auxdata.c = 0.6;

% Leader's trajectory long
    % auxdata.vl = @(t) (t<=120) * 30 ...
    %             + (((t>120) + (t<= 240)) ==2) .* (-t/6 + 50) ...
    %             + (((t>240) + (t<= 420)) ==2) .* 10 ...
    %             + (((t>420) + (t<= 540)) ==2) .* (t/6 - 60) ...
    %             + (t>540) .* 30;
auxdata.utime = (0:1:240)';
auxdata.time = (0:0.1:240)';
% Leader's trajectory short
vl = @(t) (t<=80) .* 30 ...
            + (((t>80) + (t<= 120)) ==2) .* (-t./8 + 40) ...
            + (((t>120) + (t<= 150)) ==2) .* 25 ...
            + (((t>150) + (t<= 190)) ==2) .* (t./8 + 25/4) ...
            + (t>190) .* 30;
vl = vl(auxdata.time);
vl = smoothdata(vl, "movmean", 200);
auxdata.vl = griddedInterpolant(auxdata.time, vl);

x0 = (eq_headway(auxdata.vl(0), auxdata)+auxdata.l) * flip(0:1:auxdata.len_platoon)';
auxdata.v0 = ones(auxdata.len_platoon, 1) * auxdata.vl(0); 
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[~, xl] = ode45(@(t,x) auxdata.vl(t), auxdata.time, x0(1), opts);
auxdata.x0 = x0(2:end);
auxdata.xl = griddedInterpolant(auxdata.time, xl);

%Initial solution  
U0 = diff(auxdata.vl(auxdata.utime) - 5);
U0 = [U0;0];
U0 = [U0, U0, U0]; 
[X0, V0, A0] = system_solve(U0, auxdata);

%% Run Optimizer (first iteration)

% First iteration

<<<<<<< HEAD
options = optimoptions('fmincon','Display','iter-detailed', ...
                        'SpecifyObjectiveGradient', true ,...
                        'FunValCheck','on', 'DerivativeCheck', 'off',...
                        'maxfunevals',1e6, 'StepTolerance',1e-12, ...
                        'algorithm', 'sqp', ...
                        'MaxIterations', 400);

fun = @(U) objective_gradient_acc(U, auxdata); 
nonlcon = []; %@(U) nlc(U, auxdata);
=======
%options = optimoptions('fmincon','Display','iter-detailed', ...
                        %'SpecifyObjectiveGradient',params("use_gradient"),...
                        %'FunValCheck','on', 'DerivativeCheck', 'off',...
                        %'maxfunevals',1e6, 'StepTolerance',1e-12, ...
                        %'algorithm', 'interior-point');

fun = @(U) objective_gradient(U, params, scenario); 
% Nonlinear constraints: accepts a vector or array x and returns two arrays, c(x) and ceq(x)
nonlcon = @(U) nlc(U, params, scenario);
>>>>>>> 49b6f44... ipopt implementation attempt #1
A = [];
b = [];
% [A, b] = lc(auxdata); 
Aeq = []; beq = []; 
a_min = -3 * ones(size(U0)); 
<<<<<<< HEAD
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
=======
a_max = 1.5 * ones(size(U0));

funcs.objective = @(U) U;
funcs.gradient = fun;
funcs.constraints = @(nonlcol) nonlcol;
% Jacobian
%funcs.jacobian  = @(nonlcol) jacobian(nonlcol, U);
%funcs.jacobian = @train_NLP_constraints_jacobian;
%funcs.jacobianstructure = @train_NLP_constraints_jacobian_pattern;
%funcs.jacobianstructure = @(nonlcol) jacobian(nonlcol, U);
%starting point
x0 = X0;


% Lower/Upper bounds on constraints (don't know if these are correct
option.cl = a_min;    
option.cu = a_max;

% Set the IPOPT options -These might have to change
option.ipopt.print_level           = 3;
option.ipopt.jac_c_constant        = 'yes';
option.ipopt.hessian_approximation = 'limited-memory';
option.ipopt.mu_strategy           = 'adaptive';
option.ipopt.tol                   = 1e-7;

[x info] = ipopt(x0,funcs,option);
%[U_star, f_val, ~, output, ~, grad] = fmincon(fun, U0, A, b, Aeq, beq, a_min, a_max, nonlcon, options);

% [U_star,~,~,~,~,grad,~] = fmincon(fun, U0, A,b,Aeq,beq,lb,ub,nonlcon,options);
%info;

% these need to be uncommented and added in later
%Fu = griddedInterpolant(scenario("time"),U_star);
%U_star = @(t) Fu(t);
%[X_star, V_star] = system_solve(U_star, params, scenario);
>>>>>>> 49b6f44... ipopt implementation attempt #1

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