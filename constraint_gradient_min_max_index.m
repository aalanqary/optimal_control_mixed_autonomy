function [c, dc] = constraint_gradient_min_max_index(X, V, A, Fx, Fv, Fa, Fu, auxdata, i, U_vec)
    
    % Objective 
%     [X, V, A] = system_solve(U_vec, auxdata);
    cmin = Cmin(auxdata.Ia(i), X, V, A, auxdata);
    cmax = Cmax(auxdata.Ia(i), X, V, A, auxdata);
    c = [cmin, cmax];

    % Gradient
     if nargout >= 2
         %PQ = ode3(@(t,PQ) F_adjoint_mod(P_dot, Q_dot, i, auxdata), flip(auxdata.time), PQ0);
         PQ0 = get_adjoint_ic(X, V, U_vec, auxdata);
         PQ = ode3(@(t,PQ) F_adjoint_min(t, i, PQ, Fx(t)', Fv(t)', Fu(t)',Fa(t)', auxdata), flip(auxdata.time), PQ0);
         PQ = flip(PQ,1);
         Q = PQ(:, auxdata.len_platoon+1:end);
         %P = PQ(:, 1:auxdata.len_platoon);
         Q_short = griddedInterpolant(auxdata.time, Q);
         Q_short = Q_short(auxdata.utime);
         
         % Lu - zeta fu
         dc_min = Q_short(:, auxdata.Ia) + phi_partial_min(flip(auxdata.time), 0, [], [], U_vec, [],  "u", auxdata)'; 
         
         PQ = ode3(@(t,PQ) F_adjoint_max(t, i, PQ, Fx(t)', Fv(t)', Fu(t)',Fa(t)', auxdata), flip(auxdata.time), PQ0);
         PQ = flip(PQ,1);
         Q = PQ(:, auxdata.len_platoon+1:end);
         %P = PQ(:, 1:auxdata.len_platoon);
         Q_short = griddedInterpolant(auxdata.time, Q);
         Q_short = Q_short(auxdata.utime);
         
         % Lu - zeta fu
         dc_max = Q_short(:, auxdata.Ia) + phi_partial_max(flip(auxdata.time), 0, [], [], U_vec, [],  "u", auxdata)'; 
         min_reshape = reshape(dc_min.', [], 1); % instead of m columns for num of AVs, it makes it to one column, stacking them up
         max_reshape = reshape(dc_max.', [], 1);
         dc = [min_reshape, max_reshape];
     end 
 
end 

%%%%% constraint Function %%%%%
function cmin = Cmin(i, X, V, A, auxdata)
    eps = auxdata.eps; 
    Xl = [auxdata.xl(auxdata.time), X]; 
    h = (Xl(:, i) - X(:, i) - auxdata.l) - auxdata.d_min; 
    cond1 = (h< -eps);
    cond2 = and((h<= eps), (h >= - eps));
    c = h.*cond1 + (-1/(4*eps))* (h - eps).^2 .* cond2;
    cmin = -auxdata.gamma - trapz(auxdata.time, c); 
end 

function cmax = Cmax(i, X, V, A, auxdata)
    eps = auxdata.eps; 
    Xl = [auxdata.xl(auxdata.time), X]; 
    h =  auxdata.d_max - (Xl(:, i) - X(:, i) - auxdata.l);
    cond1 = (h< -eps);
    cond2 = and((h<= eps), (h >= - eps));
    c = h.*cond1 + (-1/(4*eps))* (h - eps).^2 .* cond2;
    cmax = -auxdata.gamma - trapz(auxdata.time, c); 
end 


%%%%% Adjoint system %%%%%
function PQ_dot = F_adjoint_min(t, i, PQ, X, V, U, A, auxdata)
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
    P_dot(auxdata.Ih) = phi_partial_min(t, i, X, V, U, A, "xh", auxdata) ...
        - Q(auxdata.Ih) .* ACC_partial(Xl(auxdata.Ih), X(auxdata.Ih), Vl(auxdata.Ih), V(auxdata.Ih), "x", auxdata) ...
        - Q(auxdata.Ih+1) .* ismembc(auxdata.Ih+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ih), Xf(auxdata.Ih), V(auxdata.Ih), Vf(auxdata.Ih), "xl", auxdata);
    
    Q_dot(auxdata.Ih) = phi_partial_min(t, i, X, V, U, A, "vh", auxdata) ...
        - P(auxdata.Ih) ...
        - Q(auxdata.Ih) .* ACC_partial(Xl(auxdata.Ih), X(auxdata.Ih), Vl(auxdata.Ih), V(auxdata.Ih), "v", auxdata)  ...
        - Q(auxdata.Ih+1) .* ismembc(auxdata.Ih+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ih), Xf(auxdata.Ih), V(auxdata.Ih), Vf(auxdata.Ih), "vl", auxdata);
    
    %%%% Automated vehicles (AV) %%%%
    P_dot(auxdata.Ia) = phi_partial_min(t, i, X, V, U, A, "xa", auxdata) ...
        - Q(auxdata.Ia+1) .* ismembc(auxdata.Ia+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ia), Xf(auxdata.Ia), V(auxdata.Ia), Vf(auxdata.Ia), "xl", auxdata);
    
    Q_dot(auxdata.Ia) = phi_partial_min(t, i, X, V, U, A, "va", auxdata) ...
        - P(auxdata.Ia) ...
        - Q(auxdata.Ia+1) .* ismembc(auxdata.Ia+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ia), Xf(auxdata.Ia), V(auxdata.Ia), Vf(auxdata.Ia), "vl", auxdata);
    PQ_dot = [P_dot; Q_dot];
end 


%%%%% Partial derivatives of the running cost %%%%%
function dl = phi_partial_min(t, i, X, V, U, A, var, auxdata)
    Ia = auxdata.Ia;
    Ih = auxdata.Ih;
    if var ~= "u"
        Xl = [auxdata.xl(t); X]; 
        Xf = [X(2:end); 0]; 
        Vl = [auxdata.vl(t); V]; 
        Vf = [V(2:end); 0]; 
        Af = [A(2:end); 0]; 
    end
  
    if var == "xh" 
        car_index = auxdata.Ia(i); % check this line
        dl = 0 * Ih; 
        if ismembc(car_index-1, Ih)
            eps = auxdata.eps; 
            dh = Xl(car_index) - X(car_index) - auxdata.l - auxdata.d_min - auxdata.eps;
            dh = 2*dh*1;
            cond1 = (dh< -eps);
            cond2 = and((dh<= eps), (dh >= - eps));
            phi = ones(length(dh)).*cond1 + (-1/(4*eps))*dh.* cond2; % length of dh is 1
            idx = find(auxdata.Ih == car_index-1);
            dl(idx) = -phi; %%d phi/ dxi-1; 
        end
    elseif var == "xa"
        car_index = auxdata.Ia(i); % check this line
        dl = 0 * Ia; 
        eps = auxdata.eps;
        dh = Xl(car_index) - X(car_index) - auxdata.l - auxdata.d_min - auxdata.eps;
        dh = 2*dh*-1;
        cond1 = (dh< -eps);
        cond2 = and((dh<= eps), (dh >= - eps));
        phi = -ones(length(dh)).*cond1 + (-1/(4*eps))*dh.* cond2;
        idx = find(auxdata.Ia == car_index);
        dl(idx) = -phi; %%d phi/ dxi; 
        if ismembc(car_index-1, Ia)
            eps = auxdata.eps; 
            dh = Xl(car_index) - X(car_index) - auxdata.l - auxdata.d_min - auxdata.eps;
            dh = 2*dh*1;
            cond1 = (dh< -eps);
            cond2 = and((dh<= eps), (dh >= - eps));
            phi = ones(length(dh)).*cond1 + (-1/(4*eps))*dh.* cond2; % flips cond1
            idx = find(auxdata.Ia == car_index - 1);
            dl(idx) = -phi; %%d phi/ dxi; 
        end
    elseif var == "vh"
        dl = 0 * Ih; 
    elseif var == "va"
        dl = 0 * Ia;
    else 
       dl = 0*U;    
    end
    dl = dl';
end 


%%%%% Adjoint system %%%%%
function PQ_dot = F_adjoint_max(t, i, PQ, X, V, U, A, auxdata)
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
    P_dot(auxdata.Ih) = phi_partial_max(t, i, X, V, U, A, "xh", auxdata) ...
        - Q(auxdata.Ih) .* ACC_partial(Xl(auxdata.Ih), X(auxdata.Ih), Vl(auxdata.Ih), V(auxdata.Ih), "x", auxdata) ...
        - Q(auxdata.Ih+1) .* ismembc(auxdata.Ih+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ih), Xf(auxdata.Ih), V(auxdata.Ih), Vf(auxdata.Ih), "xl", auxdata);
    
    Q_dot(auxdata.Ih) = phi_partial_max(t, i, X, V, U, A, "vh", auxdata) ...
        - P(auxdata.Ih) ...
        - Q(auxdata.Ih) .* ACC_partial(Xl(auxdata.Ih), X(auxdata.Ih), Vl(auxdata.Ih), V(auxdata.Ih), "v", auxdata)  ...
        - Q(auxdata.Ih+1) .* ismembc(auxdata.Ih+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ih), Xf(auxdata.Ih), V(auxdata.Ih), Vf(auxdata.Ih), "vl", auxdata);
    
    %%%% Automated vehicles (AV) %%%%
    P_dot(auxdata.Ia) = phi_partial_max(t, i, X, V, U, A, "xa", auxdata) ...
        - Q(auxdata.Ia+1) .* ismembc(auxdata.Ia+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ia), Xf(auxdata.Ia), V(auxdata.Ia), Vf(auxdata.Ia), "xl", auxdata);
    
    Q_dot(auxdata.Ia) = phi_partial_max(t, i, X, V, U, A, "va", auxdata) ...
        - P(auxdata.Ia) ...
        - Q(auxdata.Ia+1) .* ismembc(auxdata.Ia+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ia), Xf(auxdata.Ia), V(auxdata.Ia), Vf(auxdata.Ia), "vl", auxdata);
    PQ_dot = [P_dot; Q_dot];
end 

%%%%% Partial derivatives of the running cost for dmax %%%%%
function dl = phi_partial_max(t, i, X, V, U, A, var, auxdata)
    Ia = auxdata.Ia;
    Ih = auxdata.Ih;
    if var ~= "u"
        Xl = [auxdata.xl(t); X]; 
        Xf = [X(2:end); 0]; 
        Vl = [auxdata.vl(t); V]; 
        Vf = [V(2:end); 0]; 
        Af = [A(2:end); 0]; 
    end
    if var == "xh" 
        car_index = auxdata.Ia(i); % check this line
        dl = 0 * Ih; 
        if ismembc(car_index-1, Ih)
            eps = auxdata.eps; 
            dh = auxdata.d_max - Xl(car_index) + X(car_index) + auxdata.l - auxdata.eps;
            dh = 2*dh*-1;
            cond1 = (dh< -eps);
            cond2 = and((dh<= eps), (dh >= - eps));
            phi = -ones(length(dh)).*cond1 + (-1/(4*eps))*dh.* cond2; % length of dh is 1
            idx = find(auxdata.Ih == car_index-1);
            dl(idx) = -phi; %%d phi/ dxi-1; 
        end
    elseif var == "xa"
        car_index = auxdata.Ia(i); % check this line
        dl = 0 * Ia; 
        eps = auxdata.eps;
        dh = auxdata.d_max - Xl(car_index) + X(car_index) + auxdata.l - auxdata.eps;
        dh = 2*dh*1;
        cond1 = (dh< -eps);
        cond2 = and((dh<= eps), (dh >= - eps));
        phi = ones(length(dh)).*cond1 + (-1/(4*eps))*dh.* cond2;
        idx = find(auxdata.Ia == car_index-1);
        dl(idx) = -phi; %%d phi/ dxi-1; 
        if ismembc(car_index-1, Ia)
            dh = auxdata.d_max - Xl(car_index) + X(car_index) + auxdata.l - auxdata.eps;
            dh = 2*dh*-1;
            cond1 = (dh< -eps);
            cond2 = and((dh<= eps), (dh >= - eps));
            phi = -ones(length(dh)).*cond1 + (-1/(4*eps))*dh.* cond2;
            idx = find(auxdata.Ia == car_index-1);
            dl(idx) = -phi; %%d phi/ dxi-1; 
        end
    elseif var == "vh"
        dl = 0 * Ih; 
    elseif var == "va"
        dl = 0 * Ia;
    else 
       dl = 0*U;    
    end
    dl = dl';
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