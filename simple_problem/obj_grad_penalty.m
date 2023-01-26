function df = obj_grad_penalty(U, auxdata)  
    [time_v, v] = system_solve(U, auxdata);
    %[dp1_dv, dp2_dv] = dp_dv(U, v, auxdata);
    %[dp1_du, dp2_du] = dp_du(U, v, auxdata);
    v = griddedInterpolant(time_v, v, "linear");
    % Interpolating dp_du
    %dp1_du_interp = griddedInterpolant(time_v, dp1_du, "linear")
    %dp2_du_interp = griddedInterpolant(time_v, dp2_du, "linear")
    lambdaT = 0;
    %length(1 + lambdaT *(auxdata.k1 + 2*auxdata.k2*v(0)) + auxdata.gamma*(dp1_du + dp2_du))
    %length(lambdaT)
    Lambda = ode5(@(t,lambda) helper(U(t), v(t), lambda, auxdata), flip(time_v), lambdaT);
    size(Lambda)
    Lambda = flip(Lambda);
    Lambda = griddedInterpolant(time_v, Lambda, "linear");
    df = arrayfun(@(a,b) -1*integral(@(t) -1 * auxdata.gamma*(dp1_du + dp2_du) + Lambda, a, b), auxdata.tau(1:end-1), auxdata.tau(2:end));  
    if isrow(df)
        df = df';
    end 
    df = [df; 0];
end

function f = dp1_dv(u, v, auxdata)
    h = auxdata.g + auxdata.k3*v^2 - u; 
    if h < -auxdata.eps
        f = 2*auxdata.k3*v;
    elseif (h >= -auxdata.eps) && (h<= auxdata.eps)
        f = (-1/4*auxdata.eps) * (4*auxdata.k3^2 * v^3 + 4*auxdata.k3*v*(auxdata.g - u - auxdata.eps));
    else
        f = 0; 
    end
end 

function f = dp2_dv(u, v, auxdata)
    h = auxdata.g + auxdata.k3*v^2 + u; 
    if h < -auxdata.eps
        f = 2*auxdata.k3*v(t);
    elseif (h >= -auxdata.eps) && (h<= auxdata.eps)
        f = (-1/4*auxdata.eps) * (4*auxdata.k3^2 * v^3 + 4*auxdata.k3*v*(auxdata.g + u - auxdata.eps));
    else
        f = 0; 
    end
end 

function h = helper(u, v, lambda, auxdata)
    h = 1 + lambda*(auxdata.k1 + 2*auxdata.k2*v) + auxdata.gamma*(dp1_dv(u, v, auxdata) + dp2_dv(u, v, auxdata));
end