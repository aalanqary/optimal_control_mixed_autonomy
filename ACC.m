function a = ACC(xl, xf, vl, vf, auxdata)
    l = auxdata.l; 
    safe_dist = auxdata.safe_dist;
    k = auxdata.k;
    v_max = auxdata.v_max;

    head_way = xl - xf - l;
    c = tanh(l + safe_dist);
    V = v_max .* ((tanh(k*head_way - safe_dist) + c)./(1 + c));
    a = auxdata.alpha *(V - vf) + ...
        auxdata.beta * (vl - vf) ./ (head_way.^2);
end 