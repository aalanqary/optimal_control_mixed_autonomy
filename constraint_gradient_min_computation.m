function [X, V, A, Fx, Fv, Fa, Fu] = constraint_gradient_min_computation(U_vec, auxdata)
    
    % Objective 
    [X, V, A] = generic_system_solve(U_vec, auxdata);
%     c = C(i, X, V, A, auxdata);
    
%     figure(1)
%     plot(U_vec)
%     title("AV Acceleration. Objective = ", z)
%     drawnow;
% 
%     figure(2)
%     [Xlin, ~] = linear_ode3_system_solve(U_vec, auxdata);
%     plot(auxdata.time, auxdata.xl(auxdata.time) - X(:, 1) - auxdata.l, auxdata.utime, auxdata.xl(auxdata.utime) - Xlin - auxdata.l)
%     title("AV Headway")
%     drawnow;
% 
%     figure(3)
%     plot([auxdata.vl(auxdata.time), V])
%     title("Velocities")
%     drawnow;


    % Gradient
    if nargout >= 2
        Fx = griddedInterpolant(auxdata.time, X);
        Fv = griddedInterpolant(auxdata.time, V);
        Fa = griddedInterpolant(auxdata.time, A);
        Fu = griddedInterpolant(auxdata.utime, U_vec, "previous");
%          [P_dot, Q_dot] = F_adjoint(t, PQ, X, V, U, A, auxdata);
%         [P_dot, Q_dot]  = F_adjoint(t, PQ, Fx(t)', Fv(t)', Fu(t)',Fa(t)', auxdata);
%         PQ = ode3(@(t,PQ) F_adjoint(t, i, PQ, Fx(t)', Fv(t)', Fu(t)',Fa(t)', auxdata), flip(auxdata.time), PQ0);
%         PQ = flip(PQ,1);
%         Q = PQ(:, auxdata.len_platoon+1:end);
%         P = PQ(:, 1:auxdata.len_platoon);
%         Q_short = griddedInterpolant(auxdata.time, Q);
%         Q_short = Q_short(auxdata.utime);
%         
%         % Lu - zeta fu
%         dc = Q_short(:, auxdata.Ia) + phi_partial(0, [], [], U_vec, [],  "u", auxdata);
%         figure(4); 
%         plot(dz)
%         drawnow
    end 
 
end 

% %%%%% constraint Function %%%%%
% function c = C(i, X, V, A, auxdata)
%     % IMPLEMENT CONSTRAINTS 
%     eps = auxdata.eps; 
%     Xl = [auxdata.xl(auxdata.time), X];
%     h = Xl(:, i) - X(:, i) - auxdata.l - auxdata.d_min; 
%     cond1 = (h< -eps);
%     cond2 = and((h<= eps), (h >= - eps));
%     c = h.*cond1 + - (1/(4*eps))* (h - eps).^2 .* cond2;
%     c = -auxdata.gamma - trapz(auxdata.time, c); 
% end 



% %%%%% Adjoint system %%%%%
 function [P_dot, Q_dot] = F_adjoint(t, PQ, X, V, U, A, auxdata)
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
     P_dot(auxdata.Ih) = - Q(auxdata.Ih) .* ACC_partial(Xl(auxdata.Ih), X(auxdata.Ih), Vl(auxdata.Ih), V(auxdata.Ih), "x", auxdata) ...
         - Q(auxdata.Ih+1) .* ismembc(auxdata.Ih+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ih), Xf(auxdata.Ih), V(auxdata.Ih), Vf(auxdata.Ih), "xl", auxdata);
     
     Q_dot(auxdata.Ih) = (0 * Ih) ...
         - P(auxdata.Ih) ...
         - Q(auxdata.Ih) .* ACC_partial(Xl(auxdata.Ih), X(auxdata.Ih), Vl(auxdata.Ih), V(auxdata.Ih), "v", auxdata)  ...
         - Q(auxdata.Ih+1) .* ismembc(auxdata.Ih+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ih), Xf(auxdata.Ih), V(auxdata.Ih), Vf(auxdata.Ih), "vl", auxdata);
     
     %%%% Automated vehicles (AV) %%%%
     P_dot(auxdata.Ia) = - Q(auxdata.Ia+1) .* ismembc(auxdata.Ia+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ia), Xf(auxdata.Ia), V(auxdata.Ia), Vf(auxdata.Ia), "xl", auxdata);
     
     Q_dot(auxdata.Ia) = (0 * Ia) ...
         - P(auxdata.Ia) ...
         - Q(auxdata.Ia+1) .* ismembc(auxdata.Ia+1, auxdata.Ih)' .* ACC_partial(X(auxdata.Ia), Xf(auxdata.Ia), V(auxdata.Ia), Vf(auxdata.Ia), "vl", auxdata);
 end 

% %%%%% Partial derivatives of the running cost %%%%%
%  function dl = phi_partial(t, i, X, V, U, A, var, auxdata)
%      Ia = auxdata.Ia;
%      Ih = auxdata.Ih;
%      if var ~= "u"
%          Xl = [auxdata.xl(t); X]; 
%          Xf = [X(2:end); 0]; 
%          Vl = [auxdata.vl(t); V]; 
%          Vf = [V(2:end); 0]; 
%          Af = [A(2:end); 0]; 
%      end
%      if var == "xh" 
%          dl = 0 * Ih; 
%          if ismembc(i-1, Ih)
%              dl(i-1) = %%d phi/ dxi-1; 
%          end
%      elseif var == "xa"
%          dl = 0 * Ia; 
%          dl(i) = %%d phi/ dxi; 
%          if ismembc(i-1, Ia)
%              dl(i-1) = %%d phi/ dxi-1; 
%          end
%      elseif var == "vh"
%          dl = 0 * Ih; 
%      elseif var == "va"
%          dl = 0 * Ia;
%      else 
%         dl = 2*U;    
%      end 
%  end 

% %%%%% Partial derivatives of the ACC function %%%%%
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