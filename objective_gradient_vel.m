function [z, dz] = objective_gradient_vel(U_vec, auxdata)
    
    % Objective 
    [X, V, A] = system_solve(U_vec, auxdata);
    z = J(X, V, A, auxdata);
    figure(1)
    plot(V)
    drawnow;
    % Gradient
    if nargout >= 2
        Fx = griddedInterpolant(auxdata.time, X);
        Fv = griddedInterpolant(auxdata.time, V);
        Fu = griddedInterpolant(auxdata.utime, U_vec, "previous");
        PQ0 = get_adjoint_ic(X, V, U_vec, auxdata);
        PQ = ode5(@(t,PQ) F_adjoint(t, PQ, Fx(t)', Fv(t)', Fu(t)', auxdata), flip(auxdata.time), PQ0);
        PQ = flip(PQ,1);
        Q = PQ(:, auxdata.len_platoon+1:end);

        X_short = Fx(auxdata.time);
        V_short = Fv(auxdata.time);
        Q_short = griddedInterpolant(auxdata.time, Q);
        Q_short = Q_short(auxdata.utime);
        
        % Lu - zeta fu
        dz = Q_short(:, auxdata.Ia) + L_partial(X_short(:, auxdata.Ia), V_short(:, auxdata.Ia), U_vec, "u", auxdata);
        figure(3); 
        plot(dz)
        drawnow
 
%         dz = Q(:, auxdata.Ia) + L_partial(X(:, auxdata.Ia), V(:, auxdata.Ia), Fu(auxdata.time), "u", auxdata);
%         dz = griddedInterpolant(auxdata.time, dz);
%         dz = arrayfun(@(a,b) integral(@(t) dz(t), a, b), auxdata.utime(1:end-1), auxdata.utime(2:end));
%         dz = [dz;0];
    end 
 
end 

%%%%% Objective Function %%%%%
function j = J(X, V, A, auxdata)
    running_cost = sum(V.^2, "all")/length(auxdata.time);
    terminal_cost = 0; %-sum(X(end, :));
    j = auxdata.mu1 * running_cost + auxdata.mu2 * terminal_cost;
end 



%%%%% Adjoint system %%%%%
function PQ_dot = F_adjoint(t, PQ, X, V, U, auxdata)
    P = PQ(1:auxdata.len_platoon); %x multipliers 
    Q = [PQ(auxdata.len_platoon+1:end); 0]; %v multipliers 
    P_dot = zeros(length(P), 1);
    Q_dot = zeros(length(Q)-1, 1);
    Xl = [auxdata.xl(t); X]; 
    Xf = [X(2:end); 0]; 
    Vl = [auxdata.vl(t); V]; 
    Vf = [V(2:end); 0]; 

    %%%% Human vehicles (HV) %%%%
    %%%% dot(zeta) = -Ly + zeta fy
    P_dot(auxdata.Ih) = -L_partial(X, V, U, "xh", auxdata) ...
        - Q(auxdata.Ih) .* ACC_partial(Xl(auxdata.Ih), X(auxdata.Ih), Vl(auxdata.Ih), V(auxdata.Ih), "x", auxdata) ...
        - Q(auxdata.Ih+1) .* ismember(auxdata.Ih+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ih), Xf(auxdata.Ih), V(auxdata.Ih), Vf(auxdata.Ih), "xl", auxdata);
    
    Q_dot(auxdata.Ih) = - L_partial(X, V, U, "vh", auxdata) ...
        - P(auxdata.Ih) ...
        - Q(auxdata.Ih) .* ACC_partial(Xl(auxdata.Ih), X(auxdata.Ih), Vl(auxdata.Ih), V(auxdata.Ih), "v", auxdata)  ...
        - Q(auxdata.Ih+1) .* ismember(auxdata.Ih+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ih), Xf(auxdata.Ih), V(auxdata.Ih), Vf(auxdata.Ih), "vl", auxdata);
    
    %%%% Automated vehicles (AV) %%%%
    P_dot(auxdata.Ia) = -L_partial(X, V, U, "xa", auxdata) ...
        - Q(auxdata.Ia+1) .* ismember(auxdata.Ia+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ia), Xf(auxdata.Ia), V(auxdata.Ia), Vf(auxdata.Ia), "xl", auxdata);
    
    Q_dot(auxdata.Ia) = -L_partial(X, V, U, "va", auxdata) ...
        - P(auxdata.Ia) ...
        - Q(auxdata.Ia+1) .* ismember(auxdata.Ia+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ia), Xf(auxdata.Ia), V(auxdata.Ia), Vf(auxdata.Ia), "vl", auxdata);
    PQ_dot = [P_dot; Q_dot];
end 

%%%%% Partial derivatives of the running cost %%%%%
function l_partial = L_partial(X, V, U, var, auxdata)
    switch var
        case "xh" 
            l_partial = zeros(size(X(auxdata.Ih)));
        case "xa"
            l_partial = zeros(size(X(auxdata.Ia)));
        case "vh"
            l_partial =  (2) * V(auxdata.Ih);
        case "va"
            l_partial =  (2) * V(auxdata.Ia); 
        case "u" 
           l_partial = zeros(size(U));            
    end 
end 

%%%%% Partial derivatives of the ACC function %%%%%
function a_partial = ACC_partial(Xl, X, Vl, V, var, auxdata)
    head_way = Xl - X - auxdata.l;
    V_prime = (auxdata.v_max * (sech(head_way - auxdata.safe_dist)).^2) ./ (1+tanh(auxdata.l+auxdata.safe_dist));
    switch var
        case "xl" 
            a_partial = auxdata.alpha * V_prime - 2 * auxdata.beta * ((Vl-V)./head_way.^3);
        case "x"
            a_partial = -1*auxdata.alpha * V_prime + 2*auxdata.beta*((Vl-V)./head_way.^3);
        case "vl"
            a_partial = auxdata.beta ./ head_way.^2;
        case "v"
            a_partial = -1*auxdata.alpha - auxdata.beta ./ head_way.^2;
    end 
end 

function PQ_0 = get_adjoint_ic(X, V, U, auxdata)
    Q_0 = zeros(auxdata.len_platoon, 1);
    P_0 = zeros(auxdata.len_platoon, 1);
    PQ_0 = [P_0; Q_0];
end 