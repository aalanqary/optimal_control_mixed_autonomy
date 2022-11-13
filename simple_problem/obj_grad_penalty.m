function df = obj_grad_penalty(U, auxdata)  
    [time_v, v] = system_solve(U, auxdata);    
    lambdaT = 0;
    v = griddedInterpolant(time_v, v, "linear");
    [dp1_dv, dp2_dv] = dp_dv(U, v, auxdata);
    [dp1_du, dp2_du] = dp_du(U, v, auxdata);
    lambdaT = 0;
    Lambda = ode5(@(t,lambda) -1 - lambda *(auxdata.k1 + 2*auxdata.k2*v(t)) - auxdata.gamma*(dp1_du + dp2_du), flip(time_v), lambdaT);    
    Lambda = flip(Lambda);
    Lambda = griddedInterpolant(time_v, Lambda, "linear");
    df = arrayfun(@(a,b) -1*integral(@(t) -1 * auxdata.gamma*(dp1_du + dp2_du) + Lambda, a, b), auxdata.tau(1:end-1), auxdata.tau(2:end));  
    if isrow(df)
        df = df';
    end 
    df = [df; 0];
end