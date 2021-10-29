function a = ACC(xl, xf, vl, vf, params)
    V = @(head_way) params("v_max") * ((tanh(head_way - params("safe_dist")) + tanh(params("l") + params("safe_dist")))./(1 + tanh(params("l") + params("safe_dist"))));
    a = params("alpha") *(V(xl - xf - params("l")) - vf) + params("beta") * (vl - vf) ./ ((xl - xf - params("l")).^2);
end 