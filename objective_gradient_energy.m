function [z, dz] = objective_gradient_energy(U_vec, auxdata)
    
    % Objective 
    [X, V, A] = system_solve(U_vec, auxdata);
    z = J(V, A, auxdata);
    
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
        V_short = griddedInterpolant(auxdata.time, V);
        V_short = V_short(auxdata.utime);
        % Lu - zeta fu
        dz = Q_short(:, auxdata.Ia) + L_partial(0, [],[],[], V_short,[],[], U_vec, [],  "u", auxdata);
    end 
 
end 

%%%%% Objective Function %%%%%
function j = J(V, A, auxdata)
    E = simplified_fuel_model(V,A,'RAV4');
    running_cost = auxdata.mu1 * sum(E.^2, "all") * auxdata.dt;    
    terminal_cost = 0; %-sum(X(end, :));
    j = running_cost + terminal_cost;
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
    P_dot(auxdata.Ih) = -L_partial(t, X, Xl, Xf, V, Vl, Vf, U, A, "xh", auxdata) ...
        - Q(auxdata.Ih) .* ACC_partial(Xl(auxdata.Ih), X(auxdata.Ih), Vl(auxdata.Ih), V(auxdata.Ih), "x", auxdata) ...
        - Q(auxdata.Ih+1) .* ismembc(auxdata.Ih+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ih), Xf(auxdata.Ih), V(auxdata.Ih), Vf(auxdata.Ih), "xl", auxdata);
    
    Q_dot(auxdata.Ih) = - L_partial(t, X, Xl, Xf, V, Vl, Vf, U, A, "vh", auxdata) ...
        - P(auxdata.Ih) ...
        - Q(auxdata.Ih) .* ACC_partial(Xl(auxdata.Ih), X(auxdata.Ih), Vl(auxdata.Ih), V(auxdata.Ih), "v", auxdata)  ...
        - Q(auxdata.Ih+1) .* ismembc(auxdata.Ih+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ih), Xf(auxdata.Ih), V(auxdata.Ih), Vf(auxdata.Ih), "vl", auxdata);
    
    %%%% Automated vehicles (AV) %%%%
    P_dot(auxdata.Ia) = -L_partial(t, X, Xl, Xf, V, Vl, Vf, U, A, "xa", auxdata) ...
        - Q(auxdata.Ia+1) .* ismembc(auxdata.Ia+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ia), Xf(auxdata.Ia), V(auxdata.Ia), Vf(auxdata.Ia), "xl", auxdata);
    
    Q_dot(auxdata.Ia) = -L_partial(t, X, Xl, Xf, V, Vl, Vf, U, A, "va", auxdata) ...
        - P(auxdata.Ia) ...
        - Q(auxdata.Ia+1) .* ismembc(auxdata.Ia+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ia), Xf(auxdata.Ia), V(auxdata.Ia), Vf(auxdata.Ia), "vl", auxdata);
    PQ_dot = [P_dot; Q_dot];
end 

%%%%% Partial derivatives of the running cost %%%%%
function dl = L_partial(t, X, Xl, Xf, V, Vl, Vf, U, A, var, auxdata)
    C = auxdata.C;
    p = auxdata.p;
    q = auxdata.q;

    Ia = auxdata.Ia;
    Ih = auxdata.Ih;
    if var ~= "u"
%         Xl = [auxdata.xl(t); X]; 
%         Xf = [X(2:end); 0]; 
%         Vl = [auxdata.vl(t); V]; 
%         Vf = [V(2:end); 0]; 
        Af = [A(2:end); 0]; 
    end
    if var == "xh" 
        da = ACC_partial(Xl(Ih), X(Ih), Vl(Ih), V(Ih), "x", auxdata);
        daf = ACC_partial(X(Ih), Xf(Ih), V(Ih), Vf(Ih), "xl", auxdata);
        dl = (p(1) + p(2).* V(Ih) + p(3) .* V(Ih).^2) .* da ...
            + (A(Ih) > 0) .* 2 .* ( q(1) +  q(2) .* V(Ih)) .* A(Ih) .* da;
        dlf = + ismembc(Ih+1, Ih)' .* (( p(1) +  p(2).* Vf(Ih) +  p(3) .* Vf(Ih).^2) .* daf ...
            + (Af(Ih) > 0) .* 2 .* ( q(1) +  q(2) * Vf(Ih)) .* Af(Ih) .* daf); 
        dl = dl + dlf;

    elseif var == "xa"
        daf = ACC_partial(X(Ia), Xf(Ia), V(Ia), Vf(Ia), "xl", auxdata);
        dl = 0;
        dlf = ismembc(Ia+1, Ih)' .* (( p(1) +  p(2).* Vf(Ia) +  p(3) .* Vf(Ia).^2) .* daf ...
            + (Af(Ia) > 0) .* 2 .* ( q(1) +  q(2) * Vf(Ia)) .* Af(Ia) .* daf); 
        dl = dl + dlf;

    elseif var == "vh"
        da = ACC_partial(Xl(Ih), X(Ih), Vl(Ih), V(Ih), "v", auxdata);
        daf = ACC_partial(X(Ih), Xf(Ih), V(Ih), Vf(Ih), "vl", auxdata);
        dl = C(2) + 2*C(3)*V(Ih) + 3*C(4)*V(Ih).^2 ...
             + ( p(1) +  p(2)*V(Ih) +  p(3)*V(Ih).^2) .* da ...
             + ( p(2) + 2* p(3)*V(Ih)) .* A(Ih) ...
             + (A(Ih) > 0) .* (2 * ( q(1) +  q(2) * V(Ih)) .* A(Ih) .* da ...
             +  q(2) * A(Ih).^2);
        dlf = + ismembc(Ih+1, Ih)' .* (( p(1) +  p(2).* Vf(Ih) +  p(3) .* Vf(Ih).^2) .* daf ...
            + (Af(Ih) > 0) .* 2 .* ( q(1) +  q(2) * Vf(Ih)) .* Af(Ih) .* daf); 
        dl = dl + dlf;
    elseif var == "va"
        daf = ACC_partial(X(Ia), Xf(Ia), V(Ia), Vf(Ia), "vl", auxdata);

        dl  = C(2) + 2*C(3)*V(Ia) + 3*C(4)*V(Ia).^2 ...
              + ( p(2) + 2* p(3)*V(Ia)) .* A(Ia);
        dlf = ismembc(Ia+1, Ih)' .* (( p(1) +  p(2).* Vf(Ia) +  p(3) .* Vf(Ia).^2) .* daf ...
            + (Af(Ia) > 0) .* 2 .* ( q(1) +  q(2) * Vf(Ia)) .* Af(Ia) .* daf); 
        dl = dl + dlf;
    else 
       dl =  p(1) +  p(2) .* V(:, Ia) +  p(3) .* V(:, Ia).^2 ...
          + (U(:, Ia) > 0) .* (2.* q(1).*U(:, Ia) + 2.* q(2).*V(:, Ia).*U(:, Ia));
    end 
end 

%%%%% Partial derivatives of the ACC function %%%%%
function a_partial = ACC_partial(Xl, X, Vl, V, var, auxdata)
    l = auxdata.l; 
    safe_dist = auxdata.safe_dist;
    k = auxdata.k;
    v_max = auxdata.v_max;

    head_way = Xl - X - l;
    V_prime = (v_max * k * (1./cosh(k*head_way - safe_dist)).^2) ./ (1+tanh(l+safe_dist));
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