function [A, b] = lc(auxdata)
    Tf = auxdata.utime(end);
    nt = length(auxdata.utime);
%     x_0 = auxdata.v0(1); 
%     X_l = auxdata.xl(auxdata.time);
    v_0 = auxdata.v0(1); 
    T = (Tf/nt) * ones(nt-1,nt);
    T = [zeros(1,nt);T];
    T = tril(T);
%     V = v_0 * ones(params("nt"), 1) + T*U;
%     X = x_0 * ones(params("nt"), 1) + v_0 * params("t_int") + T*T*U;
%     A = [T*T; -T*T;  -T]; 
%     b = [X_l - x_0 - v_0 * params("t_int") - params("l") - params("eps"); ...
%         params("gamma") - X_l + x_0 + v_0 * params("t_int") + params("l"); ...
%         v_0 * ones(params("nt"), 1)];
    A = -T;
    b = v_0 * ones(length(auxdata.utime), 1);

end