function dl = L_partial(t, X, V, U, A, var, auxdata)
    if var ~= "u"
        Xl = [auxdata.xl(t); X]; 
        Xf = [X(2:end); 0]; 
        Vl = [auxdata.vl(t); V]; 
        Vf = [V(2:end); 0]; 
        Af = [A(2:end); 0]; 
    
        h = Xl(auxdata.Ia) - X(auxdata.Ia) - auxdata.l;
        hf = [h(2:end); 0];
    end
    switch var
        case "xh" 
            dl_acc = 2 * A(auxdata.Ih) .* ACC_partial(Xl(auxdata.Ih), X(auxdata.Ih), Vl(auxdata.Ih), V(auxdata.Ih), "x", auxdata) ...
                + ismember(auxdata.Ih+1, auxdata.Ih)' * 2 .* Af(auxdata.Ih) .* ACC_partial(X(auxdata.Ih), Xf(auxdata.Ih), V(auxdata.Ih), Vf(auxdata.Ih), "xl", auxdata);
            
            dl_penalty = ismember(auxdata.Ih+1, auxdata.Ia)' * -1 .* ((auxdata.a*auxdata.b)/((auxdata.b*hf + auxdata.c).^2 + 1));

            dl = dl_acc + dl_penalty;
        case "xa"
            dl_acc = 2 * ismember(auxdata.Ia+1, auxdata.Ih)' .* Af(auxdata.Ia) .* ACC_partial(X(auxdata.Ia), Xf(auxdata.Ia), V(auxdata.Ia), Vf(auxdata.Ia), "xl", auxdata);
            
            dl_penalty = -1 .* ((auxdata.a*auxdata.b)/((auxdata.b*h + auxdata.c).^2 + 1)) .* -1 ...
                + ismember(auxdata.Ia+1, auxdata.Ia)' * -1 .* ((auxdata.a*auxdata.b)/((auxdata.b*hf + auxdata.c).^2 + 1));
          
            dl = dl_acc + dl_penalty;
        case "vh"
            dl = 2 * A(auxdata.Ih) .* ACC_partial(Xl(auxdata.Ih), X(auxdata.Ih), Vl(auxdata.Ih), V(auxdata.Ih), "v", auxdata) ...
                + 2 * ismember(auxdata.Ih+1, auxdata.Ih)' .* Af(auxdata.Ih) .* ACC_partial(X(auxdata.Ih), Xf(auxdata.Ih), V(auxdata.Ih), Vf(auxdata.Ih), "vl", auxdata);
          
        case "va"
            dl = 2 * ismember(auxdata.Ia+1, auxdata.Ih)' .* Af(auxdata.Ia) .* ACC_partial(X(auxdata.Ia), Xf(auxdata.Ia), V(auxdata.Ia), Vf(auxdata.Ia), "vl", auxdata);

        case "u" 
           dl = 2*U;         
    end 
end 