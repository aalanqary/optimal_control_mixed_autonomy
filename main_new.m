%% Define problem params
data = importfile("driving_traj/2021-03-08-22-35-14_2T3MWRFVXLW056972_masterArray_0_5373.csv");
% data = importfile("2021-04-20-21-42-34_2T3MWRFVXLW056972_masterArray_0_4395.csv");
params = containers.Map;
% ODE parameters 
time = data.Time - min(data.Time);
time = time(1:100);
params("t_0") = min(time);
params("T") = max(time); 
params("nt") = length(time);
params("t_int") = linspace(params("t_0"), params("T"), params("nt"))'; 
params("ode_opt") = [];%odeset('RelTol',1e-10,'AbsTol',1e-12);
% IDM params
params("safe_dist") = 2.5; %float
params("v_max") = 35; %float
params("l") = 4.5; %float
params("alpha") = 0.5;
params("beta") = 21;
% objective function params 
params("mu1") = 200;
params("mu2") = 1;
% optimizer params 
params("use_gradient") = true; 
params("initalize_as_IDM") = true; 

%%Define interpolation function 
global fun_interp forward_sys_solve
fun_interp = @(val, t) interp1(params("t_int"), val, t, 'nearest');%% try nearest interp or spline or cubic 

forward_sys_solve = 0;

%% Define scenario 
config = [1, 0]; %1 represents AV and 0 represents HV
vl = data.Velocity * (1000/3600); 
vl = vl(1:100);
headway_0 = 11;
x_0 = headway_0 * flip(0:length(config))'; %% Includes leader
velocity_0 = 14;
v_0 = velocity_0*ones(length(config),1);
scenario = generate_scenario(config, x_0, v_0, @(t) fun_interp(vl, t), params);
xl = scenario("x_leader");
xl = xl(params("t_int"));
scenario("number") = "5373_mu1_200_mu2_1_10";


%% Get initial guess for optimizer 
% Set initial controller guess
if params("initalize_as_IDM")
    [X_0, V_0, U_0] = initial_solution(containers.Map(keys(scenario), values(scenario)), containers.Map(keys(params), values(params)));
else
    U_0 = @(t) [0*t, 0*t, 0*t];
    U_0 = U_0(params("t_int"));
end 
% [X_1, V_1] = system_solve(U_0, params, scenario); 
% xl = scenario("x_leader");
%     xl = xl(params("t_int"));
%     vl = scenario("v_leader");
%     vl = vl(params("t_int"));
%     X_l = [xl, X_1]; 
%     V_l = [vl, V_1]; 
%     U_1 = ACC(X_l(:, 1), X_1(:, 1), ... 
%               V_l(:, 1), V_1(:, 1), params);
%% Run optimizer 
tic
options = optimoptions('fmincon','Display','iter-detailed', ...
                        'SpecifyObjectiveGradient',params("use_gradient"),...
                        'FunValCheck','on', 'DerivativeCheck', 'off',...
                        'maxfunevals',1e6, 'StepTolerance',1e-12,'algorithm', 'sqp');

fun = @(U) objective_gradient(U, params, scenario); 
nonlcon = @(U) nlc(U, params, scenario);
A1 = []; b = []; Aeq = []; beq = []; 
lb = -3 * ones(size(U_0)); 
ub = 3 * ones(size(U_0));
[U_star,~,~,~,~,grad,~] = fmincon(fun, U_0, A1,b,Aeq,beq,lb,ub,nonlcon,options);
toc
[X_star, V_star] = system_solve(U_star, params, scenario);

%%
sprintf('scenarios/scenario%s.mat', scenario("number"));
save(sprintf('scenarios/scenario%s.mat', scenario("number")), "U_star", "X_star", "V_star", "params", "scenario")

%%

%Plot initial solution 
f_scenario = figure();
plot(params("t_int"), xl, 'LineWidth',2)
hold on 
plot(params("t_int"), X_0(:, [1,2]), 'LineWidth',2)
xlabel("time (s)", 'FontSize', 18)
ylabel("Position (m)", 'FontSize', 18)
saveas(f_scenario, sprintf('scenarios/%s_pos.png', scenario("number")));

f_scenario = figure();
plot(params("t_int"), vl, 'LineWidth',2)
hold on 
plot(params("t_int"), V_0(:, [1,2]), 'LineWidth',2)
xlabel("time (s)", 'FontSize', 18)
ylabel("Velocity (m/s)", 'FontSize', 18)
saveas(f_scenario, sprintf('scenarios/%s_vel.png', scenario("number")));


%% Plot final solution 
f_scenario = figure();
plot(params("t_int"), xl, 'LineWidth',2)
hold on 
plot(params("t_int"), X_star(:, [1,2]), 'LineWidth',2)
xlabel("time (s)", 'FontSize', 18)
ylabel("Position (m)", 'FontSize', 18)
saveas(f_scenario, sprintf('scenarios/%s_posSTAR.png', scenario("number")));

f_scenario = figure();
plot(params("t_int"), vl, 'LineWidth',2)
hold on 
plot(params("t_int"), V_star(:, [1,2]), 'LineWidth',2)
xlabel("time (s)", 'FontSize', 18)
ylabel("Velocity (m/s)", 'FontSize', 18)
saveas(f_scenario, sprintf('scenarios/%s_velSTAR.png', scenario("number")));

%% Prepare results 
mkdir('scenarios', scenario("number"))
save(sprintf('scenarios/%s/results.mat', scenario("number")), "U_star", "X_star", "V_star", "params", "scenario")
fc_0 = params("T")/params("nt") * sum(FC(X_0, V_0, U_0, scenario, params), "all");
fc_star = params("T")/params("nt") * sum(FC(X_star, V_star, U_star, scenario, params), "all");
fig_0 = plot_x_v(X_0, V_0, params, scenario);
saveas(fig_0, sprintf('scenarios/%s/initial.png', scenario("number")))
fig_star = plot_x_v(X_star, V_star, params, scenario);
saveas(fig_star, sprintf('scenarios/%s/final.png', scenario("number")))

reduction = (fc_0 - fc_star)/fc_0;
