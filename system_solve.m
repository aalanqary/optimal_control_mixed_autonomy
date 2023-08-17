function [X, V, A] = system_solve(U_vec, auxdata, leader) %% TODO what is causing the long running time
    if isempty(U_vec)
        U = @(t) [];
    else
        U = griddedInterpolant(auxdata.utime, U_vec, "previous");
    end
    

    XV = ode3(@(t,XV) F(t, XV, U(t)', auxdata, leader), ... 
                     auxdata.time, [auxdata.x0; auxdata.v0]);
    V = XV(:, end - auxdata.len_platoon + 1: end);
    X = XV(:, 1:end - auxdata.len_platoon);
    

    if nargout > 2
        A = zeros(size(X));
        Xl = [leader.x(auxdata.time), X]; 
        Vl = [leader.v(auxdata.time), V];  
        Ah = ACC(Xl(:, auxdata.Ih), X(:, auxdata.Ih), ... 
                 Vl(:, auxdata.Ih), V(:, auxdata.Ih), auxdata); 
        A(:, auxdata.Ih) = Ah;
        if ~isempty(U)
            A(:, auxdata.Ia) = U(auxdata.time); 
        end 
    end 
%     Linear solver
%     x_0 = scenario("x_0");
%     x_0 = x_0(1); 
%     v_0 = scenario("v_0");
%     v_0 = v_0(1); 
%     T = (params("T")/params("nt")) * ones(params("nt")-1,params("nt"));
%     T = [zeros(1,params("nt"));T];
%     T = tril(T);
%     V = v_0 * ones(params("nt"), 1) + T*U;
%     X = x_0 * ones(params("nt"), 1) + v_0 * params("t_int") + T*T*U;

end 

function XV_dot = F(t, XV, U, auxdata, leader) 
    X = XV(1:auxdata.len_platoon);
    V = XV(auxdata.len_platoon+1:end);
    Xl = [leader.x(t); X]; 
    Vl = [leader.v(t); V]; 
    a = ACC(Xl(auxdata.Ih), X(auxdata.Ih), ...
            Vl(auxdata.Ih), V(auxdata.Ih), auxdata);
    X_dot = V;
    V_dot = zeros(length(V), 1);
    V_dot(auxdata.Ih) = a;
    V_dot(auxdata.Ia) = U; 
    XV_dot = [X_dot; V_dot];
end 