function [z, dz] = objective_gradient_acc_penalty(U_vec, auxdata)
    
    % Objective 
    [X, V, A] = generic_system_solve(U_vec, auxdata);
    z = J(X, V, A, auxdata);
    
    % Gradient
    if nargout >= 2
        Fx = griddedInterpolant(auxdata.time, X);
        Fv = griddedInterpolant(auxdata.time, V);
        Fa = griddedInterpolant(auxdata.time, A);
        Fu = griddedInterpolant(auxdata.utime, U_vec, "previous");
        PQ0 = get_adjoint_ic(X, V, U_vec, auxdata);
        PQ = ode3(@(t,PQ) F_adjoint(t, PQ, Fx(t)', Fv(t)', Fu(t)',Fa(t)', auxdata), flip(auxdata.time), PQ0);
        PQ = flip(PQ,1);
        Q = PQ(:, auxdata.len_platoon+1:end);
        P = PQ(:, 1:auxdata.len_platoon);
        Q_short = griddedInterpolant(auxdata.time, Q);
        Q_short = Q_short(auxdata.utime);
        
        % Lu - zeta fu
        dz = Q_short(:, auxdata.Ia) + L_partial(0, [], [], U_vec, [],  "u", auxdata);
%         figure(4); 
%         plot(dz)
%         drawnow
    end 
 
end 

%%%%% Objective Function %%%%%
function j = J(X, V, A, auxdata)
    % Compute headways of AVs
    Xl = [auxdata.xl(auxdata.time), X]; 
    Vl = [auxdata.vl(auxdata.time), V]; 
    h = Xl(:, auxdata.Ia) - X(:, auxdata.Ia) - auxdata.l;
    
    % arctan new barrier function (a(-arctan(bx+c)+pi/2) might have to use
    % sym(pi)
%     auxdata.mu2_tmp = 1./(1+exp(-auxdata.iter+2));
    running_cost = auxdata.mu1 * 0.1 * sum(A.^2, "all") ...
        + auxdata.mu2 * 0.1 * sum(auxdata.a.*(-atan(auxdata.b.*(h) + auxdata.c) + pi/2), "all") ...
        + auxdata.mu3 * 0.1 * sum(auxdata.d.*(atan(auxdata.e.*(h) + auxdata.f) + pi/2), "all"); 
    terminal_cost = 0; %-sum(X(end, :));
    j = running_cost + terminal_cost;
    %       figure(3)  
    %       plot(Vl)
    %       legend('v1', 'v2', 'v3', 'v4', 'v5');
    %     title("Velocity")

%     drawnow;
%     display(auxdata.iter)
%     if mod(auxdata.iter,5) == 0
%          figure(3)
%          plot(Vl)
%          title("Velocity")
%          drawnow;
%     end
%     auxdata.iter = auxdata.iter + 1;

    %       figure(6)
    %       plot(A)
    %       title("Acceleration, Objective = ", j)
    %     drawnow;

%     penalty2 = auxdata.mu2 .* auxdata.a.*(-atan(auxdata.b.*h + auxdata.c) + pi/2);
%     penalty3 = auxdata.mu3 * auxdata.d.*(atan(auxdata.e.*h + auxdata.f) + pi/2);
    %       figure(4)
    %       plot(Xl(:, 1) - Xl(:, 2) - auxdata.l)
    %       title("AV Headway leader first AV")
    %       figure(5)
    %       plot(Xl(:, 3) - Xl(:, 4) - auxdata.l)
    %       title("AV Headway second AV")
%     figure(4)
%     plot(penalty2)
%     title("Penalty 2")
%     legend('AV1', 'AV2')
%     figure(14)
%     plot(penalty3)
%     title("Penalty 3")
%     legend('AV1', 'AV2')
%     drawnow;
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
    Ia = auxdata.Ia;
    Ih = auxdata.Ih;
    if var ~= "u"
        Xl = [auxdata.xl(t); X]; 
        Xf = [X(2:end); 0]; 
        Vl = [auxdata.vl(t); V]; 
        Vf = [V(2:end); 0]; 
        Af = [A(2:end); 0]; 
    
        h = Xl(1:end-1) - X - auxdata.l;
        hf = [h(2:end); 0];
    end
        if var == "xh" 
            dl_acc = 2 * A(Ih) .* ACC_partial(Xl(Ih), X(Ih), Vl(Ih), V(Ih), "x", auxdata) ...
            + ismembc(Ih+1, Ih)' * 2 .* Af(Ih) .* ACC_partial(X(Ih), Xf(Ih), V(Ih), Vf(Ih), "xl", auxdata);            
            
            dl_min = ismembc(auxdata.Ih+1, auxdata.Ia)' * -1 .* ((auxdata.a*auxdata.b)./((auxdata.b*hf(auxdata.Ih) + auxdata.c).^2 + 1));
            dl_max = ismembc(auxdata.Ih+1, auxdata.Ia)' * 1 .* ((auxdata.d*auxdata.e)./((auxdata.e*hf(auxdata.Ih) + auxdata.f).^2 + 1));
            
            dl = auxdata.mu1.*dl_acc + auxdata.mu2.*dl_min + auxdata.mu3*dl_max;
        elseif var == "xa"
            dl_acc = 2 * ismembc(Ia+1, Ih)' .* Af(Ia) .* ACC_partial(X(Ia), Xf(Ia), V(Ia), Vf(Ia), "xl", auxdata);            
            dl_min = -1 .* ((auxdata.a*auxdata.b)./((auxdata.b*h(Ia) + auxdata.c).^2 + 1)) .* -1 ...
                + ismembc(Ia+1, Ia)' * -1 .* ((auxdata.a*auxdata.b)./((auxdata.b*hf(Ia) + auxdata.c).^2 + 1));
            
            dl_max = 1 .* ((auxdata.d*auxdata.e)./((auxdata.e*h(Ia) + auxdata.f).^2 + 1)) .* -1 ...
                + ismember(Ia+1, Ia)' * 1 .* ((auxdata.d*auxdata.e)./((auxdata.e*hf(Ia) + auxdata.f).^2 + 1));
            dl = auxdata.mu1.*dl_acc + auxdata.mu2.*dl_min + auxdata.mu3*dl_max;
        elseif var == "vh"
            dl = 2 * A(Ih) .* ACC_partial(Xl(Ih), X(Ih), Vl(Ih), V(Ih), "v", auxdata) ...
            + 2 * ismembc(Ih+1, Ih)' .* Af(Ih) .* ACC_partial(X(Ih), Xf(Ih), V(Ih), Vf(Ih), "vl", auxdata);
          
        elseif var == "va"
            dl = 2 * ismembc(Ia+1, Ih)' .* Af(Ia) .* ACC_partial(X(Ia), Xf(Ia), V(Ia), Vf(Ia), "vl", auxdata);

        else 
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