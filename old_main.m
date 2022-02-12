%%Define problem params 
params = containers.Map;
% ODE parameters 
params("t_0") = 0;
params("T") = 70; 
params("nt") = 500;
params("t_int") = linspace(params("t_0"), params("T"), params("nt"))'; 
params("ode_opt") = []; %odeset('RelTol',1e-8,'AbsTol',1e-10);
% IDM params
params("safe_dist") = 2.5; %float
params("v_max") = 35; %float
params("l") = 4.5; %float
params("alpha") = 0.5;
params("beta") = 21;
% objective function params 
params("des_v") = @(t) 0*t + 21;
% optimizer params 
params("use_gradient") = true; 
params("initalize_as_IDM") = true; 
% Define interpolation function 
global fun_interp
fun_interp = @(val, t) interp1(params("t_int"), val, t);

%%Define scenario 
v_leader = @(t) (t<=20) .* (0.*t + 21) ...
                + ((20<t) & (t<=30)).* (-(3/5) .* t + 33) ...
                + ((30<t) & (t<=40)) .* (0.*t + 15) ...
                + ((40<t) & (t<=50)) .* ((3/5)*t - 9) ...
                + (50 < t) .* (0.*t + 21); 

% v_leader = @(t) (t<=10) .* (0.*t + 21) + ((10<t) & (t<=12)).* (-1 .* t +31) + ((12<t) & (t<=15)) .* (0.*t + 19) + ((15<t) & (t<=17)) .* (t + 4) + (17 < t) .* (0.*t + 21); 

            %                 + ((20<t) & (t<=35)).* (-(11/15) .* t +107/3) ...
            %                 + ((65<t) & (t<=80)) .* ((11/15)*t - (113/3)) ...

% v_leader = @(t) 0*t + 21;
config = [1, 0, 0, 0, 0]; %1 represent AV and 0 HV
headway_0 = 10;
v_0 = 21;
t_int = linspace(params("t_0"), params("T"), 2*params("nt"));
[~, x_leader] = ode45(@(t, dd) v_leader(t), t_int, length(config)*headway_0);
x_leader = @(t) interp1(t_int, x_leader, t);
scenario = generate_scenario(config, headway_0, v_0, x_leader, v_leader);
scenario("number") = "energy2_grad_nt200";

%%Plot leader velocity  
f_scenario = figure('Color','white', 'Unit', 'inches', 'Position', [0,0,5,5]);
plot(params("t_int"), v_leader(params("t_int")), 'LineWidth',2)
xlabel("time (s)", 'FontSize', 18)
ylabel("Leader velocity", 'FontSize', 18)
% saveas(f_scenario, sprintf('scenarios/scenario%s.png', scenario("number")));

%% Get initial guess for optimizer 
% Set initial controller guess
if params("initalize_as_IDM")
    [~, ~, U_0] = initial_solution(containers.Map(keys(scenario), values(scenario)), containers.Map(keys(params), values(params)));
else
    U_0 = @(t) 0*t + 21;
    U_0 = U_0(params("t_int"));
end 
% Solve the system with initial guess
[X_0, V_0] = system_solve(U_0, params, scenario);
% Plot velocity and position 
% plot_x_v(X_0, V_0, U_0, params, scenario, "initial_sol");

%% Run optimizer 
tic
options = optimoptions('fmincon','Display','iter-detailed', ...
                        'SpecifyObjectiveGradient',params("use_gradient"),...
                        'FunValCheck','on', 'DerivativeCheck', 'off',...
                        'maxfunevals',1e4, 'StepTolerance',1e-12,'algorithm', 'sqp');

fun = @(U) objective_gradient(U, params, scenario);


A1 = [];
b = [];
Aeq = [];
beq = [];
lb = zeros(size(U_0)); 
ub = [];
nonlcon = [];
[U_star,~,~,~,~,grad,~] = fmincon(fun, U_0, A1,b,Aeq,beq,lb,ub,nonlcon,options);
toc
%%

[X_star, V_star] = system_solve(U_star, params, scenario);
% Plot velocity and position 
plot_x_v(X_star, V_star, U_star, params, scenario, "final_sol");
% 
% 
% sprintf('scenarios/scenario%s.mat', scenario("number"))
% save(sprintf('scenarios/scenario%s.mat', scenario("number")), "U_star", "X_star", "V_star", "params", "scenario")
