function [X0, V0, U0] = initial_solution(auxdata)
      original_Ia = auxdata.Ia;
      original_Ih = auxdata.Ih;
      original_platoon = auxdata.platoon;
      original_time = auxdata.time;

      auxdata.platoon = zeros(1, auxdata.len_platoon);
      auxdata.Ia = find(auxdata.platoon);
      auxdata.Ih = find(auxdata.platoon - 1);
      auxdata.l = 2 * auxdata.l;
      auxdata.time = auxdata.utime; 

      [X0, V0, U0] = system_solve([], auxdata);
      auxdata.platoon = original_platoon;
      auxdata.Ia = original_Ia;
      auxdata.Ih = original_Ih;
      auxdata.l = auxdata.l/10;
      auxdata.time = original_time;
      U0 = U0(:, original_Ia);
end 