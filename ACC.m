function a = ACC(xl, xf, vl, vf, params)
    V = @(head_way) params("v_max") * ((tanh(head_way - params("safe_dist")) + tanh(params("l") + params("safe_dist")))./(1 + tanh(params("l") + params("safe_dist"))));
    head_way = xl - xf - params("l");
    a = params("alpha") *(V(head_way) - vf) + params("beta") * (vl - vf) ./ (head_way.^2);
end 