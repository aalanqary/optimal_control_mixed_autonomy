function h = eq_headway(v, auxdata)
    C = tanh(auxdata.l + auxdata.safe_dist);
    h = atanh((v * (1+C)./auxdata.v_max) - C) + auxdata.safe_dist;
    h = h/auxdata.k;
end 