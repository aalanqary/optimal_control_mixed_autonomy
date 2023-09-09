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
init = 1; 
schedule = 1; 
if schedule == 1
    mu_min = 1; 
    mu_max = 1; 
    factor = 10; 
end 


%% get initial guess
platoon = [zeros(1,20)];
results_in = "results/real_traj/init" + init + "/";
save_res = true;
[auxdata, leader] = problem_auxdata(platoon, const, traj);
U_initial = diff(smoothdata(leader.v(auxdata.utime), 'movmean', 3)) ./ diff(auxdata.utime);
U_initial = [U_initial;0] + 0.00001;
display(U_initial)
fun = @(U) objective_gradient_accel_penalty(U, auxdata, leader); 
nonlcon = []; 
% [A, b] = lc_v(auxdata); 
A = []; 
b = []; 
Aeq = []; beq = []; 
a_min = [];
a_max = []; 
tic
%[U_1av, ~, ~, ~, ~, ~] = fmincon(fun, U_initial, A, b, Aeq, beq, a_min, a_max, nonlcon, options);
[X1av, V1av, A1av] = system_solve([], auxdata, leader);
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
end 

data.objective_val = trapz(auxdata.time, sum(A1av.^2, 2));
data.optim_system_energy = sum(trapz(auxdata.time, simplified_fuel_model(V1av, A1av, 'RAV4'))); % why is this having the sum after an integration
data.min_headway = min(min(headway - auxdata.d_min, 0), [], "all");
data.max_headway = max(max(headway - auxdata.d_max, 0), [], "all");
display("energy consumption = " +  data.optim_system_energy)
display("optimization time = " + timee)
display("Minimum headway violation = " + data.min_headway)
display("Maximum headway violation = " + data.max_headway)
display("objective value = " + data.objective_val)

save(results_in + 'AV0' + penalty_iter + '.mat', 'data')


%% Solve platoon
init = 2;
platoon = [1,zeros(1,3),1,zeros(1, 3),1,zeros(1,3),1,zeros(1,3),1,zeros(1,3)]; % next is 7
[auxdata, leader] = problem_auxdata(platoon, const, traj);
platoon_name = length(auxdata.Ia) + "av_" + length(auxdata.Ih)/length(auxdata.Ia) + "hv";
results_in = "results/real_traj/init3" + "/" + platoon_name + "_" + schedule +"/"; 
save_res = true;
func_eval_count = 0;
tic
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
        elseif init == 2
            load("results/real_traj/init3/4av_4hv_1/U_10.mat")
            load("results/real_traj/init1/U_1av.mat")
            U0 = [U_star, U_star(:, 4)];
        elseif init == 3
            load("results/real_traj/init1/"+platoon_name+"_1/U_9.mat")
            U0 = U_star;
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
    [X_star, V_star, A_star] = system_solve(U_star, auxdata, leader);
    
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
    end 

    data.min_violation = min(violations);
    data.max_violation = max(violations);
    data.objective_val = trapz(auxdata.time, sum(A_star.^2, 2));
    data.optim_system_energy = sum(trapz(auxdata.time, simplified_fuel_model(V_star, A_star, 'RAV4'))); % why is this having the sum after an integration
    func_eval_count = func_eval_count + output.funcCount;
    data.func_count = func_eval_count;
    save(results_in + 'Data_' + penalty_iter + '.mat', 'data')

    display("Minimum violation = " + data.min_violation)
    display("Maximum violation = " + data.max_violation)
    display("objective value = " + data.objective_val)
    display("energy = " + data.optim_system_energy)
    auxdata.mu_min = auxdata.mu_min * factor; 
    auxdata.mu_max = auxdata.mu_max * factor; 
end 
timee = toc;
fin_val.timee = timee;
fin_val.func_count = func_eval_count;
display("total optimization time: " + fin_val.timee)
display("function counter: " + fin_val.func_count)
save(results_in + 'Total' + '.mat', 'fin_val')


