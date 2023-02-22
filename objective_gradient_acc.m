function [z, dz] = objective_gradient_acc(U_vec, auxdata)
    
    % Objective 
    [X, V, A] = system_solve(U_vec, auxdata);
    z = J(X, V, A, auxdata);
    
    % Gradient
    if nargout >= 2
        Fx = griddedInterpolant(auxdata.time, X);
        Fv = griddedInterpolant(auxdata.time, V);
        Fa = griddedInterpolant(auxdata.time, A);
        Fu = griddedInterpolant(auxdata.utime, U_vec, "previous");
        PQ0 = get_adjoint_ic(X, V, U_vec, auxdata);
        PQ = ode5(@(t,PQ) F_adjoint(t, PQ, Fx(t)', Fv(t)', Fu(t)',Fa(t)', auxdata), flip(auxdata.time), PQ0);
        PQ = flip(PQ,1);
        Q = PQ(:, auxdata.len_platoon+1:end);

        X_short = Fx(auxdata.time);
        V_short = Fv(auxdata.time);
        Q_short = griddedInterpolant(auxdata.time, Q);
        Q_short = Q_short(auxdata.utime);
        
        % Lu - zeta fu
        dz = Q_short(:, auxdata.Ia) + L_partial(0, [], [], U_vec, [],  "u", auxdata);
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
    % Compute headways of AVs
    Xl = [auxdata.xl(auxdata.time), X]; 
    Vl = [auxdata.vl(auxdata.time), V]; 
    h = Xl(:, auxdata.Ia) - X(:, auxdata.Ia) - auxdata.l;
    
    running_cost = auxdata.mu1 * sum(A.^2, "all")/length(auxdata.time) ...
        + auxdata.mu2 * sum(log(max(min(1, h), 1e-10)).^2, "all")/length(auxdata.time); %use smooth function 
    terminal_cost = 0; %-sum(X(end, :));
    j = running_cost + terminal_cost;

    figure(1)
    plot(A)
    title("Objective = ", sum(A.^2, "all")/length(auxdata.time))
    drawnow;
    
    figure(2)
    plot(Xl)
%     title("Objective = ", sum(log(max(min(1, h), 0)).^2, "all")/length(auxdata.time), "all")
    drawnow;

end 



%%%%% Adjoint system %%%%%
function PQ_dot = F_adjoint(t, PQ, X, V, U, A, auxdata)
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
    P_dot(auxdata.Ih) = -L_partial(t, X, V, U, A, "xh", auxdata) ...
        - Q(auxdata.Ih) .* ACC_partial(Xl(auxdata.Ih), X(auxdata.Ih), Vl(auxdata.Ih), V(auxdata.Ih), "x", auxdata) ...
        - Q(auxdata.Ih+1) .* ismember(auxdata.Ih+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ih), Xf(auxdata.Ih), V(auxdata.Ih), Vf(auxdata.Ih), "xl", auxdata);
    
    Q_dot(auxdata.Ih) = - L_partial(t, X, V, U, A, "vh", auxdata) ...
        - P(auxdata.Ih) ...
        - Q(auxdata.Ih) .* ACC_partial(Xl(auxdata.Ih), X(auxdata.Ih), Vl(auxdata.Ih), V(auxdata.Ih), "v", auxdata)  ...
        - Q(auxdata.Ih+1) .* ismember(auxdata.Ih+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ih), Xf(auxdata.Ih), V(auxdata.Ih), Vf(auxdata.Ih), "vl", auxdata);
    
    %%%% Automated vehicles (AV) %%%%
    P_dot(auxdata.Ia) = -L_partial(t, X, V, U, A, "xa", auxdata) ...
        - Q(auxdata.Ia+1) .* ismember(auxdata.Ia+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ia), Xf(auxdata.Ia), V(auxdata.Ia), Vf(auxdata.Ia), "xl", auxdata);
    
    Q_dot(auxdata.Ia) = -L_partial(t, X, V, U, A, "va", auxdata) ...
        - P(auxdata.Ia) ...
        - Q(auxdata.Ia+1) .* ismember(auxdata.Ia+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ia), Xf(auxdata.Ia), V(auxdata.Ia), Vf(auxdata.Ia), "vl", auxdata);
    PQ_dot = [P_dot; Q_dot];
end 

%%%%% Partial derivatives of the running cost %%%%%
function dl = L_partial(t, X, V, U, A, var, auxdata)
    if var ~= "u"
        Xl = [auxdata.xl(t); X]; 
        Xf = [X(2:end); 0]; 
        Vl = [auxdata.vl(t); V]; 
        Vf = [V(2:end); 0]; 
        Af = [A(2:end); 0]; 
    
        h = Xl(auxdata.Ia) - X(auxdata.Ia) - auxdata.l;
        hf = [h(2:end); 0];
    end
    switch var
        case "xh" 
            dl_acc = 2 * A(auxdata.Ih) .* ACC_partial(Xl(auxdata.Ih), X(auxdata.Ih), Vl(auxdata.Ih), V(auxdata.Ih), "x", auxdata) ...
                + ismember(auxdata.Ih+1, auxdata.Ih)' * 2 .* Af(auxdata.Ih) .* ACC_partial(X(auxdata.Ih), Xf(auxdata.Ih), V(auxdata.Ih), Vf(auxdata.Ih), "xl", auxdata);
            
            dl_penality = ismember(auxdata.Ih+1, auxdata.Ia)' * 2 .* (log(max(min(1, hf), 1e-10))./(max(min(1, hf), 1e-10))) .* (((hf < 1) + (hf>= 1e-10)) ==2); 
            dl = dl_acc + dl_penality;
        case "xa"
            dl_acc = 2 * ismember(auxdata.Ia+1, auxdata.Ih)' .* Af(auxdata.Ia) .* ACC_partial(X(auxdata.Ia), Xf(auxdata.Ia), V(auxdata.Ia), Vf(auxdata.Ia), "xl", auxdata);
            
            dl_penality = 2 * (log(max(min(1, h), 1e-10))./(max(min(1, h), 1e-10))) .* -1 * ((h < 1) && (h>= 1e-10)) ...
                + ismember(auxdata.Ia+1, auxdata.Ia)' * 2 .* (log(max(min(1, hf), 1e-10))./(max(min(1, hf), 1e-10))) .* (((hf < 1) + (hf>= 1e-10)) ==2);
            dl = dl_acc + dl_penality;
        case "vh"
            dl = 2 * A(auxdata.Ih) .* ACC_partial(Xl(auxdata.Ih), X(auxdata.Ih), Vl(auxdata.Ih), V(auxdata.Ih), "v", auxdata) ...
                + 2 * ismember(auxdata.Ih+1, auxdata.Ih)' .* Af(auxdata.Ih) .* ACC_partial(X(auxdata.Ih), Xf(auxdata.Ih), V(auxdata.Ih), Vf(auxdata.Ih), "vl", auxdata);
          
        case "va"
            dl = 2 * ismember(auxdata.Ia+1, auxdata.Ih)' .* Af(auxdata.Ia) .* ACC_partial(X(auxdata.Ia), Xf(auxdata.Ia), V(auxdata.Ia), Vf(auxdata.Ia), "vl", auxdata);

        case "u" 
           dl = 2*U;         
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