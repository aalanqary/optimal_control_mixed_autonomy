% penalty_iter = 1;
options = optimoptions('fmincon','Display','iter-detailed', ...
                        'SpecifyObjectiveGradient', true ,...
                        'SpecifyConstraintGradient', false ,...
                        'FunValCheck','on', 'DerivativeCheck', 'off',...
                        'maxfunevals',1e6, 'StepTolerance',1e-12, ...
                        'algorithm', 'sqp', ...
                        'MaxIterations', 200);
traj = "real";
const = "penalty_minmax";
experiment = "sequential";
init = 1; 
schedule = 1; 
if schedule == 1
    mu_min = 1; 
    mu_max = 1; 
    factor = 10; 
end 


%% get initial guess
platoon = [1, zeros(1,19)];
results_in = "results/real_traj/" + experiment + "/";
if not(isfolder(results_in))
        mkdir(results_in)
end
save_res = true;
[auxdata, leader] = problem_auxdata(platoon, const, traj);
U_initial = diff(smoothdata(leader.v(auxdata.utime), 'movmean', 3)) ./ diff(auxdata.utime);
U_initial = [U_initial;0] + 0.00001;
fun = @(U) objective_gradient_accel_penalty(U, auxdata, leader); 
nonlcon = []; 
% [A, b] = lc_v(auxdata); 
A = []; 
b = []; 
Aeq = []; beq = []; 
a_min = [];
a_max = []; 
tic
[U_1av, ~, ~, ~, ~, ~] = fmincon(fun, U_initial, A, b, Aeq, beq, a_min, a_max, nonlcon, options);
[X_star, V_star, A_star] = system_solve(U_1av, auxdata, leader);
timee = toc; 

headway = leader.x(auxdata.time) - X1av - auxdata.l;
if save_res
    figure()
    plot(auxdata.time, leader.v(auxdata.time))
    hold on 
    plot(auxdata.time, V1av)
    ylabel("Velocity")
    savefig(results_in + 'velocity_1av.fig')
    
    figure()
    plot(auxdata.time, leader.x(auxdata.time))
    hold on 
    plot(auxdata.time, X1av)
    ylabel("Position")
    savefig(results_in + 'position_1av.fig')
    
    figure()
    plot(auxdata.time, headway)
    ylabel("Headway")
    savefig(results_in + 'headway_1av.fig')

    save(results_in + 'U_1av.mat', 'U_1av')
    save(results_in + 'auxdata_1av.mat', 'auxdata')
    % Save all results
    save(results_in + "all_results.mat")
end 

% data.timee = timee;
% data.min_violation = min(violations);
% data.max_violation = max(violations);
% data.objective_val = trapz(auxdata.time, sum(A_star.^2, 2));
% data.optim_system_energy = sum(trapz(auxdata.time, simplified_fuel_model(V_star, A_star, 'RAV4'))); % why is this having the sum after an integration
% save(results_in + 'Data_' + penalty_iter + '.mat', 'data')

display("optimization time = " + timee)
display("Minimum headway violation = " + min(min(headway - auxdata.d_min, 0), [], "all"))
display("Maximum headway violation = " + max(max(headway - auxdata.d_max, 0), [], "all"))
display("objective value = " + trapz(auxdata.time, sum(A1av.^2, 2)))

%% Solve platoon
platoon = [1,zeros(1, 3)];
av_num = 1;
if av_num == 1
    [auxdata, leader] = problem_auxdata(platoon, const, traj, "");
end
results_in = "results/real_traj/" + experiment + "/" + "AV" + av_num + "/";
if not(isfolder(results_in))
        mkdir(results_in)
end 
save_res = true;

for penalty_iter = 1:1:10
    if const == "penalty_minmax"
        fun = @(U) objective_gradient_accel_penalty(U, auxdata, leader); 
    elseif const == "smooth_penalty"
        fun = @(U) objective_gradient_accel_smooth_penalty(U, auxdata, leader); 
    end 
%     options = optimoptions(options, 'MaxIterations', 10 + 5*penalty_iter);
    if penalty_iter == 1
        if init == 1
            load("results/real_traj/init"+init+"/U_1av.mat")
            U0 = repmat(U_1av, [1, length(auxdata.Ia)]);
%         elseif init == 2
%             load("results/real_traj/init1/3av_3hv_1/U_10.mat")
%             load("results/real_traj/init1/U_1av.mat")
%             U0 = [U_star, U_1av];
%         elseif init == 3
%             load("results/real_traj/init1/"+platoon_name+"_1/U_9.mat")
%             U0 = U_star;
        end 
        [X0, V0, A0] = system_solve(U0, auxdata, leader); 
        Xl = [leader.x(auxdata.time), X0]; 
        headway = Xl(:,auxdata.Ia) - X0(:, auxdata.Ia) - auxdata.l;
        min_headway_violations = min(headway - auxdata.d_min, 0);
        min_headway_violations = nonzeros(min_headway_violations);
        max_headway_violations = max(headway - auxdata.d_max, 0);
        max_headway_violations = nonzeros(max_headway_violations);
        violations = [min_headway_violations; max_headway_violations];
        display("Initial objective value = " + trapz(auxdata.time, sum(A0.^2, 2)))
        display("Initial minimum violation = " + min(violations))
        display("Initial maximum violation = " + max(violations))
        auxdata.mu_min = mu_min; 
        auxdata.mu_max = mu_max; 
    else
        U0 = U_star;    
    end 
    nonlcon = []; 
    A = []; b = []; Aeq = []; beq = []; 
    a_min = []; a_max = []; 
    tic
    [U_star, f_val, ~, output, ~, grad] = fmincon(fun, U0, A, b, Aeq, beq, a_min, a_max, nonlcon, options);
    [X_star, V_star, A_star] = system_solve(U_star, auxdata, leader);
    timee = toc; 
    
    Xl = [leader.x(auxdata.time), X_star]; 
    headway = Xl(:,auxdata.Ia) - X_star(:, auxdata.Ia) - auxdata.l;
    min_headway_violations = min(headway - auxdata.d_min, 0);
    min_headway_violations = nonzeros(min_headway_violations);
    max_headway_violations = max(headway - auxdata.d_max, 0);
    max_headway_violations = nonzeros(max_headway_violations);
    violations = [min_headway_violations; max_headway_violations];
    if save_res
        figure()
        plot(auxdata.time, leader.v(auxdata.time), "color", "black")
        hold on 
        plot(auxdata.time, V_star(:, auxdata.Ih), "color", "blue", "LineWidth", 0.5)
        hold on 
        plot(auxdata.time, V_star(:, auxdata.Ia), "color", "red", "LineWidth", 1.5)
        ylabel("Velocity")
        savefig(results_in + 'velocity_' + penalty_iter + '.fig')
        
        figure()
        plot(auxdata.time, leader.x(auxdata.time), "color", "black")
        hold on 
        plot(auxdata.time, X_star(:, auxdata.Ih), "color", "blue", "LineWidth", 0.5)
        hold on 
        plot(auxdata.time, X_star(:, auxdata.Ia), "color", "red","LineWidth", 1.5)
        ylabel("Position")
        savefig(results_in + 'position_' + penalty_iter + '.fig')
        
        figure()
        plot(auxdata.time, headway)
        legend("AV ID = " + string(auxdata.Ia))
        ylabel("AV headway")
        savefig(results_in + 'headway_' + penalty_iter + '.fig')
        
    %     figure()
    %     histogram(violations, [min(violations):1:0, 0:1:max(violations)])
    %     title("Headway violations")
    %     savefig(results_in + 'violations_' + penalty_iter + '.fig')
    
        if penalty_iter == 1
            save(results_in + 'leader.mat', 'leader')
        end 
        save(results_in + 'auxadata_' + penalty_iter + '.mat', 'auxdata')
        save(results_in + 'U_' + penalty_iter + '.mat', 'U_star')
        % Save all results
        save(results_in + "all_results.mat")
    end 
    
    data.timee = timee;
    data.min_violation = min(violations);
    data.max_violation = max(violations);
    data.objective_val = trapz(auxdata.time, sum(A_star.^2, 2));
    data.optim_system_energy = sum(trapz(auxdata.time, simplified_fuel_model(V_star, A_star, 'RAV4'))); % why is this having the sum after an integration
    save(results_in + 'Data_' + penalty_iter + '.mat', 'data')

    display("optimization time = " + timee)
    display("Minimum violation = " + min(violations))
    display("Maximum violation = " + max(violations))
    display("objective value = " + trapz(auxdata.time, sum(A_star.^2, 2)))
    auxdata.mu_min = auxdata.mu_min * factor; 
    auxdata.mu_max = auxdata.mu_max * factor; 
end 

%% Sequential
% penalty_iter = 1;
clearvars;
options = optimoptions('fmincon','Display','iter-detailed', ...
                        'SpecifyObjectiveGradient', true ,...
                        'SpecifyConstraintGradient', false ,...
                        'FunValCheck','on', 'DerivativeCheck', 'off',...
                        'maxfunevals',1e6, 'StepTolerance',1e-12, ...
                        'algorithm', 'sqp', ...
                        'MaxIterations', 200);
traj = "real";
const = "penalty_minmax";
experiment = "sequential_new";
init = 1; 
schedule = 1; 
if schedule == 1
    mu_min = 1; 
    mu_max = 1; 
    factor = 10; 
end 

platoon = [1,zeros(1, 3)];

av_num = 3;

if av_num == 1
    [auxdata, leader] = problem_auxdata(platoon, const, traj, "", true, av_num);
else
    [auxdata, leader] = problem_auxdata(platoon, const, traj, "results/real_traj/" + experiment + "/" + "AV" + (av_num-1) + "/" + "all_results.mat", true, av_num);
end

results_in = "results/real_traj/" + experiment + "/" + "AV" + av_num + "/";
if not(isfolder(results_in))
        mkdir(results_in)
end 

save_res = true;

for penalty_iter = 1:1:10
    if const == "penalty_minmax"
        fun = @(U) objective_gradient_accel_penalty(U, auxdata, leader); 
    elseif const == "smooth_penalty"
        fun = @(U) objective_gradient_accel_smooth_penalty(U, auxdata, leader); 
    end 
%     options = optimoptions(options, 'MaxIterations', 10 + 5*penalty_iter);
    if penalty_iter == 1
        if init == 1
            %load("results/real_traj/init"+init+"/U_1av.mat")
            %U0 = repmat(U_1av, [1, length(auxdata.Ia)]);
            U_initial = diff(smoothdata(leader.v(auxdata.utime), 'movmean', 3)) ./ diff(auxdata.utime);
            U_initial = [U_initial;0] + 0.00001;
            U0 = U_initial;
        end 
        [X0, V0, A0] = system_solve(U0, auxdata, leader); 
        Xl = [leader.x(auxdata.time), X0]; 
        headway = Xl(:,auxdata.Ia) - X0(:, auxdata.Ia) - auxdata.l;
        min_headway_violations = min(headway - auxdata.d_min, 0);
        min_headway_violations = nonzeros(min_headway_violations);
        max_headway_violations = max(headway - auxdata.d_max, 0);
        max_headway_violations = nonzeros(max_headway_violations);
        violations = [min_headway_violations; max_headway_violations];
        display("Initial objective value = " + trapz(auxdata.time, sum(A0.^2, 2)))
        display("Initial minimum violation = " + min(violations))
        display("Initial maximum violation = " + max(violations))
        auxdata.mu_min = mu_min; 
        auxdata.mu_max = mu_max; 
    else
        U0 = U_star;    
    end 
    nonlcon = []; 
    A = []; b = []; Aeq = []; beq = []; 
    a_min = []; a_max = []; 
    tic
    [U_star, f_val, ~, output, ~, grad] = fmincon(fun, U0, A, b, Aeq, beq, a_min, a_max, nonlcon, options);
    [X_star, V_star, A_star] = system_solve(U_star, auxdata, leader);
    timee = toc; 
    
    Xl = [leader.x(auxdata.time), X_star]; 
    headway = Xl(:,auxdata.Ia) - X_star(:, auxdata.Ia) - auxdata.l;
    min_headway_violations = min(headway - auxdata.d_min, 0);
    min_headway_violations = nonzeros(min_headway_violations);
    max_headway_violations = max(headway - auxdata.d_max, 0);
    max_headway_violations = nonzeros(max_headway_violations);
    violations = [min_headway_violations; max_headway_violations];
    if save_res
        figure()
        plot(auxdata.time, leader.v(auxdata.time), "color", "black")
        hold on 
        plot(auxdata.time, V_star(:, auxdata.Ih), "color", "blue", "LineWidth", 0.5)
        hold on 
        plot(auxdata.time, V_star(:, auxdata.Ia), "color", "red", "LineWidth", 1.5)
        ylabel("Velocity")
        savefig(results_in + 'velocity_' + penalty_iter + '.fig')
        close()
        
        figure()
        plot(auxdata.time, leader.x(auxdata.time), "color", "black")
        hold on 
        plot(auxdata.time, X_star(:, auxdata.Ih), "color", "blue", "LineWidth", 0.5)
        hold on 
        plot(auxdata.time, X_star(:, auxdata.Ia), "color", "red","LineWidth", 1.5)
        ylabel("Position")
        savefig(results_in + 'position_' + penalty_iter + '.fig')
        close()

        figure()
        %plot(auxdata.time, leader.x(auxdata.time), "color", "black")
        %hold on 
        plot(auxdata.time, A_star(:, auxdata.Ih), "color", "blue", "LineWidth", 0.5)
        hold on 
        plot(auxdata.time, A_star(:, auxdata.Ia), "color", "red","LineWidth", 1.5)
        ylabel("Acceleration(Vehicles Only)")
        savefig(results_in + 'acceleration_' + penalty_iter + '.fig')
        close()
        
        figure()
        plot(auxdata.time, headway)
        legend("AV ID = " + string(auxdata.Ia))
        ylabel("AV headway")
        savefig(results_in + 'headway_' + penalty_iter + '.fig')
        close()
        
    %     figure()
    %     histogram(violations, [min(violations):1:0, 0:1:max(violations)])
    %     title("Headway violations")
    %     savefig(results_in + 'violations_' + penalty_iter + '.fig')
    
        if penalty_iter == 1
            save(results_in + 'leader.mat', 'leader')
        end 
        save(results_in + 'auxadata_' + penalty_iter + '.mat', 'auxdata')
        save(results_in + 'U_' + penalty_iter + '.mat', 'U_star')
        % Save all results
        save(results_in + "all_results.mat")
    end 
    
    data.timee = timee;
    data.min_violation = min(violations);
    data.max_violation = max(violations);
    data.objective_val = trapz(auxdata.time, sum(A_star.^2, 2));
    data.optim_system_energy = sum(trapz(auxdata.time, simplified_fuel_model(V_star, A_star, 'RAV4'))); % why is this having the sum after an integration
    save(results_in + 'Data_' + penalty_iter + '.mat', 'data')

    display("optimization time = " + timee)
    display("Minimum violation = " + min(violations))
    display("Maximum violation = " + max(violations))
    display("objective value = " + trapz(auxdata.time, sum(A_star.^2, 2)))
    auxdata.mu_min = auxdata.mu_min * factor; 
    auxdata.mu_max = auxdata.mu_max * factor; 
end 


%% Sequential platoons

clearvars;
options = optimoptions('fmincon','Display','iter-detailed', ...
                        'SpecifyObjectiveGradient', true ,...
                        'SpecifyConstraintGradient', false ,...
                        'FunValCheck','on', 'DerivativeCheck', 'off',...
                        'maxfunevals',1e6, 'StepTolerance',1e-12, ...
                        'algorithm', 'sqp', ...
                        'MaxIterations', 200);
traj = "real";
const = "penalty_minmax";
experiment = "sequential_platoons_new";
init = 1; 
schedule = 1; 
if schedule == 1
    mu_min = 1; 
    mu_max = 1; 
    factor = 10; 
end 

% CHANGE
platoon = [1,zeros(1, 3), 1,zeros(1, 3), 1,zeros(1, 3), 1,zeros(1, 3), 1,zeros(1, 3)];

% CHANGE
av_num = 5;

[auxdata, leader] = problem_auxdata(platoon, const, traj, "", true, av_num, true);

results_in = "results/real_traj/" + experiment + "/" + "AV" + av_num + "/";
if not(isfolder(results_in))
        mkdir(results_in)
end 

save_res = true;

%NEED TO UPDATE ON EVERY RUN
U0 = [];
for c = 1:av_num
    load("results/real_traj/sequential_new/" + "AV" + c + "/U_10.mat");
    U0 = [U0, U_star];
end

[X_star, V_star, A_star] = system_solve(U0, auxdata, leader); 

Xl = [leader.x(auxdata.time), X_star]; 
headway = Xl(:,auxdata.Ia) - X_star(:, auxdata.Ia) - auxdata.l;
min_headway_violations = min(headway - auxdata.d_min, 0);
min_headway_violations = nonzeros(min_headway_violations);
max_headway_violations = max(headway - auxdata.d_max, 0);
max_headway_violations = nonzeros(max_headway_violations);
violations = [min_headway_violations; max_headway_violations];

if save_res
    figure()
    plot(auxdata.time, leader.v(auxdata.time), "color", "black")
    hold on 
    plot(auxdata.time, V_star(:, auxdata.Ih), "color", "blue", "LineWidth", 0.5)
    hold on 
    plot(auxdata.time, V_star(:, auxdata.Ia), "color", "red", "LineWidth", 1.5)
    ylabel("Velocity")
    savefig(results_in + 'velocity.fig')
    close()
    
    figure()
    plot(auxdata.time, leader.x(auxdata.time), "color", "black")
    hold on 
    plot(auxdata.time, X_star(:, auxdata.Ih), "color", "blue", "LineWidth", 0.5)
    hold on 
    plot(auxdata.time, X_star(:, auxdata.Ia), "color", "red","LineWidth", 1.5)
    ylabel("Position")
    savefig(results_in + 'position.fig')
    close()
    
    figure()
    plot(auxdata.time, headway)
    legend("AV ID = " + string(auxdata.Ia))
    ylabel("AV headway")
    savefig(results_in + 'headway.fig')
    close()

    save(results_in + 'auxdata.mat', 'auxdata')
    save(results_in + 'U.mat', 'U_star')
    % Save all results
    save(results_in + "all_results.mat")
end

data.min_violation = min(violations);
data.max_violation = max(violations);
data.objective_val = trapz(auxdata.time, sum(A_star.^2, 2));
data.optim_system_energy = sum(trapz(auxdata.time, simplified_fuel_model(V_star, A_star, 'RAV4'))); % why is this having the sum after an integration
save(results_in + 'Data.mat', 'data')

%display("optimization time = " + timee)
display("Minimum violation = " + min(violations))
display("Maximum violation = " + max(violations))
display("Objective value = " + data.objective_val)
display("System Energy = " + data.optim_system_energy)
auxdata.mu_min = auxdata.mu_min * factor; 
auxdata.mu_max = auxdata.mu_max * factor; 

%% Coordinated Platoons with Initializations
clearvars;
options = optimoptions('fmincon','Display','iter-detailed', ...
                        'SpecifyObjectiveGradient', true ,...
                        'SpecifyConstraintGradient', false ,...
                        'FunValCheck','on', 'DerivativeCheck', 'off',...
                        'maxfunevals',1e6, 'StepTolerance',1e-12, ...
                        'algorithm', 'sqp', ...
                        'MaxIterations', 200);
traj = "real";
const = "penalty_minmax";
experiment = "coordinated";
init = 1; 
schedule = 1; 
if schedule == 1
    %CHANGED
    mu_min = 0.001; 
    %CHANGED
    mu_max = 0.001;
    %CHANGED
    factor = sqrt(10); 
end 

%CHANGE
platoon = [1,zeros(1, 3), 1,zeros(1, 3), 1,zeros(1, 3), 1,zeros(1, 3), 1,zeros(1, 3)];

%CHANGE
av_num = 5;

[auxdata, leader] = problem_auxdata(platoon, const, traj, "");

results_in = "results/real_traj/" + experiment + "/" + "AV" + av_num + "/";
if not(isfolder(results_in))
        mkdir(results_in)
end 

save_res = true;

% Initialize time and function eval metrics
tic
func_evals = 0;

for penalty_iter = 1:1:10
    if const == "penalty_minmax"
        fun = @(U) objective_gradient_accel_penalty(U, auxdata, leader); 
    elseif const == "smooth_penalty"
        fun = @(U) objective_gradient_accel_smooth_penalty(U, auxdata, leader); 
    end 
%     options = optimoptions(options, 'MaxIterations', 10 + 5*penalty_iter);
    if penalty_iter == 1
        U0 = [];
        for c = 1:av_num
            load("results/real_traj/sequential/" + "AV" + c + "/U_10.mat");
            U0 = [U0, U_star];
        end 
        [X0, V0, A0] = system_solve(U0, auxdata, leader); 
        Xl = [leader.x(auxdata.time), X0]; 
        headway = Xl(:,auxdata.Ia) - X0(:, auxdata.Ia) - auxdata.l;
        min_headway_violations = min(headway - auxdata.d_min, 0);
        min_headway_violations = nonzeros(min_headway_violations);
        max_headway_violations = max(headway - auxdata.d_max, 0);
        max_headway_violations = nonzeros(max_headway_violations);
        violations = [min_headway_violations; max_headway_violations];
        display("Initial objective value = " + trapz(auxdata.time, sum(A0.^2, 2)))
        display("Initial minimum violation = " + min(violations))
        display("Initial maximum violation = " + max(violations))
        auxdata.mu_min = mu_min; 
        auxdata.mu_max = mu_max; 
    else
        U0 = U_star;    
    end 
    nonlcon = []; 
    A = []; b = []; Aeq = []; beq = []; 
    a_min = []; a_max = []; 
    [U_star, f_val, ~, output, ~, grad] = fmincon(fun, U0, A, b, Aeq, beq, a_min, a_max, nonlcon, options);
    % Update function evaluations
    func_evals = func_evals + output.funcCount;
    [X_star, V_star, A_star] = system_solve(U_star, auxdata, leader);
    timee = toc; 
    
    Xl = [leader.x(auxdata.time), X_star]; 
    headway = Xl(:,auxdata.Ia) - X_star(:, auxdata.Ia) - auxdata.l;
    min_headway_violations = min(headway - auxdata.d_min, 0);
    min_headway_violations = nonzeros(min_headway_violations);
    max_headway_violations = max(headway - auxdata.d_max, 0);
    max_headway_violations = nonzeros(max_headway_violations);
    violations = [min_headway_violations; max_headway_violations];
    if save_res
        figure()
        plot(auxdata.time, leader.v(auxdata.time), "color", "black")
        hold on 
        plot(auxdata.time, V_star(:, auxdata.Ih), "color", "blue", "LineWidth", 0.5)
        hold on 
        plot(auxdata.time, V_star(:, auxdata.Ia), "color", "red", "LineWidth", 1.5)
        ylabel("Velocity")
        savefig(results_in + 'velocity_' + penalty_iter + '.fig')
        close()
        
        figure()
        plot(auxdata.time, leader.x(auxdata.time), "color", "black")
        hold on 
        plot(auxdata.time, X_star(:, auxdata.Ih), "color", "blue", "LineWidth", 0.5)
        hold on 
        plot(auxdata.time, X_star(:, auxdata.Ia), "color", "red","LineWidth", 1.5)
        ylabel("Position")
        savefig(results_in + 'position_' + penalty_iter + '.fig')
        close()

        figure()
        %plot(auxdata.time, leader.x(auxdata.time), "color", "black")
        %hold on 
        plot(auxdata.time, A_star(:, auxdata.Ih), "color", "blue", "LineWidth", 0.5)
        hold on 
        plot(auxdata.time, A_star(:, auxdata.Ia), "color", "red","LineWidth", 1.5)
        ylabel("Acceleration(Vehicles Only)")
        savefig(results_in + 'acceleration_' + penalty_iter + '.fig')
        close()
        
        figure()
        plot(auxdata.time, headway)
        legend("AV ID = " + string(auxdata.Ia))
        ylabel("AV headway")
        savefig(results_in + 'headway_' + penalty_iter + '.fig')
        close()
        
    %     figure()
    %     histogram(violations, [min(violations):1:0, 0:1:max(violations)])
    %     title("Headway violations")
    %     savefig(results_in + 'violations_' + penalty_iter + '.fig')
    
        if penalty_iter == 1
            save(results_in + 'leader.mat', 'leader')
        end 
        save(results_in + 'auxadata_' + penalty_iter + '.mat', 'auxdata')
        save(results_in + 'U_' + penalty_iter + '.mat', 'U_star')
        % Save all results
        save(results_in + "all_results.mat")
    end 
    
    data.timee = timee;
    data.min_violation = min(violations);
    data.max_violation = max(violations);
    data.objective_val = trapz(auxdata.time, sum(A_star.^2, 2));
    data.optim_system_energy = sum(trapz(auxdata.time, simplified_fuel_model(V_star, A_star, 'RAV4'))); % why is this having the sum after an integration
    save(results_in + 'Data_' + penalty_iter + '.mat', 'data')

    display("optimization time = " + timee)
    display("Minimum violation = " + min(violations))
    display("Maximum violation = " + max(violations))
    display("objective value = " + trapz(auxdata.time, sum(A_star.^2, 2)))
    auxdata.mu_min = auxdata.mu_min * factor; 
    auxdata.mu_max = auxdata.mu_max * factor; 
end 
timee = toc;
time_metrics.time = timee;
time_metrics.func_evals = func_evals;
save(results_in + 'time_metrics.mat', 'data')



%% Sequential Stats Combined
% Combine accelerations
load("results/real_traj/coordinated/" + "AV1/" + "all_results.mat", "A_star");
A_star_1 = A_star;
load("results/real_traj/coordinated/" + "AV2/" + "all_results.mat", "A_star");
A_star_2 = A_star;
load("results/real_traj/coordinated/" + "AV3/" + "all_results.mat", "A_star");
A_star_3 = A_star;
load("results/real_traj/coordinated/" + "AV4/" + "all_results.mat", "A_star");
A_star_4 = A_star;
load("results/real_traj/coordinated/" + "AV5/" + "all_results.mat", "A_star");
A_star_5 = A_star;

A_star = cat(2, A_star_1, A_star_2);
A_star_2av = cat(2, A_star_1, A_star_2);
A_star = cat(2, A_star, A_star_3);
A_star = cat(2, A_star, A_star_4);
A_star = cat(2, A_star, A_star_5);
save("results/real_traj/coordinated/" + "combined_stats_A_star", 'A_star');

% Combine velocities
load("results/real_traj/coordinated/" + "AV1/" + "all_results.mat", "V_star");
V_star_1 = V_star;
load("results/real_traj/coordinated/" + "AV2/" + "all_results.mat", "V_star");
V_star_2 = V_star;
load("results/real_traj/coordinated/" + "AV3/" + "all_results.mat", "V_star");
V_star_3 = V_star;
load("results/real_traj/coordinated/" + "AV4/" + "all_results.mat", "V_star");
V_star_4 = V_star;
load("results/real_traj/coordinated/" + "AV5/" + "all_results.mat", "V_star");
V_star_5 = V_star;

V_star = cat(2, V_star_1, V_star_2);
V_star_2av = cat(2, V_star_1, V_star_2);
V_star = cat(2, V_star, V_star_3);
V_star = cat(2, V_star, V_star_4);
V_star = cat(2, V_star, V_star_5);
save("results/real_traj/coordinated/" + "combined_stats_V_star", 'V_star');
%% Calculate 2 AV Combined Stats

sequential_objective_val_2av = trapz(auxdata.time, sum(A_star_2av.^2, 2))
sequential_optim_system_energy_2av = sum(trapz(auxdata.time, simplified_fuel_model(V_star_2av, A_star_2av, 'RAV4')))
%% Calculate Combined Sequential Energy Stats
load("results/real_traj/coordinated/" + "combined_stats_A_star")
load("results/real_traj/coordinated/" + "combined_stats_V_star")

sequential_objective_val = trapz(auxdata.time, sum(A_star.^2, 2))
sequential_optim_system_energy = sum(trapz(auxdata.time, simplified_fuel_model(V_star, A_star, 'RAV4')))