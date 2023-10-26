function [z, dz] = objective_gradient_accel_penalty(U_vec, auxdata, leader)
    
    % Objective 
    [X, V, A] = system_solve(U_vec, auxdata, leader);
    z = J(X, A, auxdata, leader);

    % Gradient
    if nargout >= 2
        Fx = griddedInterpolant(auxdata.time, X);
        Fv = griddedInterpolant(auxdata.time, V);
        Fa = griddedInterpolant(auxdata.time, A);
        Fu = griddedInterpolant(auxdata.utime, U_vec, "previous");
        PQ0 = get_adjoint_ic(X, V, U_vec, auxdata);
        PQ = ode3(@(t,PQ) F_adjoint(t, PQ, Fx(t)', Fv(t)', Fu(t)',Fa(t)', auxdata, leader), flip(auxdata.time), PQ0);
        PQ = flip(PQ,1);
        Q = PQ(:, auxdata.len_platoon+1:end);
%         P = PQ(:, 1:auxdata.len_platoon);

%         Q_short = griddedInterpolant(auxdata.time, Q);
%         Q_short = Q_short(auxdata.utime);
%         Q_short = Q_short(:, auxdata.Ia);

        Fz = griddedInterpolant(auxdata.time, Q(:, auxdata.Ia) + L_partial(0, [], [], Fu, [],  "u", auxdata, leader));
        Fz = @(t) Fz(t); 
        dz = zeros(size(U_vec));
        for i = 1:length(auxdata.utime)-1
            dz(i, :) = integral(Fz, auxdata.utime(i), auxdata.utime(i+1), "ArrayValued", true);
        end 

        % Lu - zeta fu
%         dz = Q_short + L_partial(0, [], [], U_vec, [],  "u", auxdata, leader);
    end 
 
end 

%%%%% Objective Function %%%%%
function j = J(X, A, auxdata, leader)
    [phi_min, phi_max] = penalty_minmax(X, auxdata, leader);
    running_cost = sum(A.^2, 2) + phi_min + phi_max;
    terminal_cost = 0; 
    j = trapz(auxdata.time, running_cost) + terminal_cost;
end 

function [phi_min, phi_max] = penalty_minmax(X, auxdata, leader)
    Xl = [leader.x(auxdata.time), X]; 
    headway = Xl(:, auxdata.Ia) - X(:, auxdata.Ia) - auxdata.l; 
    phi_min = auxdata.mu_min * sum(min(headway-auxdata.d_min, 0).^2, 2);
    phi_max = auxdata.mu_max * sum(max(headway-auxdata.d_max, 0).^2, 2);
end 




%%%%% Adjoint system %%%%%
function PQ_dot = F_adjoint(t, PQ, X, V, U, A, auxdata, leader)
    P = PQ(1:auxdata.len_platoon); %x multipliers 
    Q = [PQ(auxdata.len_platoon+1:end); 0]; %v multipliers 
    P_dot = zeros(length(P), 1);
    Q_dot = zeros(length(Q)-1, 1);
    Xl = [leader.x(t); X]; 
    Xf = [X(2:end); 0]; 
    Vl = [leader.v(t); V]; 
    Vf = [V(2:end); 0]; 

    %%%% Human vehicles (HV) %%%%
    %%%% dot(zeta) = -Ly + zeta fy
    P_dot(auxdata.Ih) = -L_partial(t, X, V, U, A, "xh", auxdata, leader) ...
        - Q(auxdata.Ih) .* ACC_partial(Xl(auxdata.Ih), X(auxdata.Ih), Vl(auxdata.Ih), V(auxdata.Ih), "x", auxdata) ...
        - Q(auxdata.Ih+1) .* ismembc(auxdata.Ih+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ih), Xf(auxdata.Ih), V(auxdata.Ih), Vf(auxdata.Ih), "xl", auxdata);
    
    Q_dot(auxdata.Ih) = - L_partial(t, X, V, U, A, "vh", auxdata, leader) ...
        - P(auxdata.Ih) ...
        - Q(auxdata.Ih) .* ACC_partial(Xl(auxdata.Ih), X(auxdata.Ih), Vl(auxdata.Ih), V(auxdata.Ih), "v", auxdata)  ...
        - Q(auxdata.Ih+1) .* ismembc(auxdata.Ih+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ih), Xf(auxdata.Ih), V(auxdata.Ih), Vf(auxdata.Ih), "vl", auxdata);
    
    %%%% Automated vehicles (AV) %%%%
    P_dot(auxdata.Ia) = -L_partial(t, X, V, U, A, "xa", auxdata, leader) ...
        - Q(auxdata.Ia+1) .* ismembc(auxdata.Ia+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ia), Xf(auxdata.Ia), V(auxdata.Ia), Vf(auxdata.Ia), "xl", auxdata);
    
    Q_dot(auxdata.Ia) = -L_partial(t, X, V, U, A, "va", auxdata, leader) ...
        - P(auxdata.Ia) ...
        - Q(auxdata.Ia+1) .* ismembc(auxdata.Ia+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ia), Xf(auxdata.Ia), V(auxdata.Ia), Vf(auxdata.Ia), "vl", auxdata);
    PQ_dot = [P_dot; Q_dot];
end 

%%%%% Partial derivatives of the running cost %%%%%
function dl = L_partial(t, X, V, U, A, var, auxdata, leader)
    Ia = auxdata.Ia;
    Ih = auxdata.Ih;
    if var ~= "u"
        Xl = [leader.x(t); X]; 
        Xf = [X(2:end); 0]; 
        Vl = [leader.v(t); V]; 
        Vf = [V(2:end); 0]; 
        Af = [A(2:end); 0]; 
        headway_min = min(Xl(1:end-1) - X - auxdata.l - auxdata.d_min, 0); 
        headway_min = [headway_min; 0];
        headway_max = max(Xl(1:end-1) - X - auxdata.l - auxdata.d_max, 0); 
        headway_max = [headway_max; 0];
    end
    if var == "xh" 
        dl_acc = 2 * A(Ih) .* ACC_partial(Xl(Ih), X(Ih), Vl(Ih), V(Ih), "x", auxdata) ...
            + ismembc(Ih+1, Ih)' * 2 .* Af(Ih) .* ACC_partial(X(Ih), Xf(Ih), V(Ih), Vf(Ih), "xl", auxdata);
        dl_penalty = 2 * ismembc(Ih+1, Ia)' .* (auxdata.mu_min * headway_min(Ih+1) + auxdata.mu_max * headway_max(Ih+1)); 
        dl = dl_acc + dl_penalty;
    elseif var == "xa"
        dl_acc = 2 * ismembc(Ia+1, Ih)' .* Af(Ia) .* ACC_partial(X(Ia), Xf(Ia), V(Ia), Vf(Ia), "xl", auxdata); 
        dl_penalty = - 2 * (auxdata.mu_min * headway_min(Ia) + auxdata.mu_max * headway_max(Ia)) ...
                     + 2 * ismembc(Ia+1, Ia)' .* (auxdata.mu_min * headway_min(Ia+1) + auxdata.mu_max * headway_max(Ia+1));
        dl = dl_acc + dl_penalty;
    elseif var == "vh"
        dl = 2 * A(Ih) .* ACC_partial(Xl(Ih), X(Ih), Vl(Ih), V(Ih), "v", auxdata) ...
            + 2 * ismembc(Ih+1, Ih)' .* Af(Ih) .* ACC_partial(X(Ih), Xf(Ih), V(Ih), Vf(Ih), "vl", auxdata);
    elseif var == "va"
        dl = 2 * ismembc(Ia+1, Ih)' .* Af(Ia) .* ACC_partial(X(Ia), Xf(Ia), V(Ia), Vf(Ia), "vl", auxdata);
    else 
       dl = 2*U(auxdata.time);    
    end 
end 

%%%%% Partial derivatives of the ACC function %%%%%
function a_partial = ACC_partial(Xl, X, Vl, V, var, auxdata)
    l = auxdata.l; 
    safe_dist = auxdata.safe_dist;
    k = auxdata.k;
    v_max = auxdata.v_max;

    head_way = Xl - X - l;
    V_prime = (v_max * k * (sech(k*head_way - safe_dist)).^2) ./ (1+tanh(l+safe_dist));
    if var == "xl" 
            a_partial = k * auxdata.alpha * V_prime - 2 * auxdata.beta * ((Vl-V)./head_way.^3);
    elseif var == "x"
            a_partial = -k * auxdata.alpha * V_prime + 2*auxdata.beta*((Vl-V)./head_way.^3);
    elseif var == "vl"
            a_partial = auxdata.beta ./ head_way.^2;
    elseif var == "v"
            a_partial = -auxdata.alpha - auxdata.beta ./ head_way.^2;
    end 
end 

function PQ_0 = get_adjoint_ic(X, V, U, auxdata)
    Q_0 = zeros(auxdata.len_platoon, 1);
    P_0 = zeros(auxdata.len_platoon, 1);
    PQ_0 = [P_0; Q_0];
end 