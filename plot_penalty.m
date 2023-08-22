function [] = plot_penalty(auxdata)
    h = linspace(0,100);
    out = sum((atan(auxdata.k.*h + auxdata.l) + pi/2), "all");
    figure(4)
    plot(h, out)
    drawnow;
end

