function a_partial = ACC_partial(xl, xf, vl, vf, var, params)
    
    V_prime = @(head_way) (params("v_max") * (sech(head_way - params("safe_dist"))).^2) ./ (1+tanh(params("l")+params("safe_dist")));
    head_way = xl - xf - params("l");
    switch var
        case 1 
            a_partial = params("alpha") * V_prime(head_way) - 2*params("beta")*((vl-vf)./head_way.^3);
        case 2
            a_partial = -1*params("alpha") * V_prime(head_way) + 2*params("beta")*((vl-vf)./head_way.^3);
        case 3
            a_partial = params("beta") ./ head_way.^2;
        case 4
            a_partial = -1*params("alpha") - params("beta") ./ head_way.^2;
    end 
end 