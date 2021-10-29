% Define scenario
params = containers.Map;
params("t_0") = 0;
params("T") = 50;
params("nt") = 500;
params("t_int") = linspace(params("t_0"), params("T"), params("nt"))'; 
params("u_0") = @(t) [9 + 0*t, 5 + 0*t];
params("safe_dist") = 2.5; %float
params("v_max") = 35; %float
params("l") = 4.5; %float
params("alpha") = 0.5;
params("beta") = 20;
% params("des_v") = @(t) (t<5) * 10 + (t>=5) * 2 ; %float
params("des_v") = @(t) 0*t + 21;

% v_leader = @(t) tan(0.1*t-1.5) + 15;
v_leader = @(t) (t<=10) .* (0.*t + 21) + ((10<t) & (t<=12)).* (-1 .* t +31) + ((12<t) & (t<=15)) .* (0.*t + 19) + ((15<t) & (t<=17)) .* (t + 4) + (17 < t) .* (0.*t + 21); 
t_int = linspace(params("t_0"), params("T"), 2*params("nt"));
[~, x_leader] = ode45(@(t, dd) v_leader(t), t_int, 6*10);
x_leader = @(t) interp1(linspace(params("t_0"), params("T"), 2*params("nt")), x_leader, t);

scenario = generate_scenario([1, 0, 0, 0, 0, 0], 10, 21, x_leader, v_leader);
scenario("number") = 2;

f_scenario = figure('Color','none', 'Unit', 'inches', 'Position', [0,0,5,5]);
plot(linspace(0, 50, 500), v_leader(linspace(0, 50, 500)), 'LineWidth',2)
xlabel("time (s)", 'FontSize', 18)
ylabel("Leader velocity", 'FontSize', 18)
saveas(f_scenario, sprintf('scenarios/scenario%i.png', scenario("number")))
%%
useGradient = true;
initalize_as_IDM = true; 

global fun_interp
fun_interp = @(val, t) interp1(params("t_int"), val, t);

% Set initial controller and visualize the system 
if initalize_as_IDM
    [X, V, U_0] = initial_solution(containers.Map(keys(scenario), values(scenario)), containers.Map(keys(params), values(params)));
else
    U_0 = params("u_0");
    U_0 = U_0(params("t_int"));
end 
system_solve(U_0, params, scenario, true, "initial_solution");

%%
options = optimoptions('fmincon','Display','iter-detailed', ...
                        'SpecifyObjectiveGradient',useGradient,...
                        'FunValCheck','on', 'DerivativeCheck', 'off',...
                        'maxfunevals',1e4, 'StepTolerance',1e-10,'algorithm', 'sqp');

fun = @(U) objective_gradient(U, scenario, params);


A = [];
b = [];
Aeq = [];
beq = [];
lb = zeros(size(U_0));
ub = [];
% ub = 40*ones(size(U_0));
nonlcon = [];
[U_star,~,~,~,~,grad,~] = fmincon(fun, U_0, A,b,Aeq,beq,lb,ub,nonlcon,options);

system_solve(U_star, params, true, "optimal_solution");