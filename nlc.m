function [c, ceq] = nlc(U_vec, auxdata)

    [X,V] = system_solve(U_vec, auxdata);
    xl = auxdata.xl(auxdata.time);
    X_l = [xl, X]; 
    headway = X_l(:, auxdata.Ia) - X(:, auxdata.Ia) - auxdata.l;
    cx = - headway; 
    cv = -V(:, auxdata.Ia);
    
%     cx = sum((min(0, headway)).^2);
%     cv = sum((min(0, V(:, auxdata.Ia))).^2);
%     
    c = [cx; cv];
    ceq = [];
end 