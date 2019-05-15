function [ts, cs] = runhm(H, c0, t, tstep)
% Run a Hamiltonian by Schrodinger's Equation:
% idc/dt = Hc;
% Similar to runh(), but this uses expm() instead of ode().
% Args:
% c0: Initial vector of c;
% t: Total time;
% tstep: Output step;
%
% Returns are similar to ode() in matlab;
% ts: List of time;
% cs: List of c at different times, with each row as c(t);
% 


    ts = 0:tstep:t;
    Ustep = expm(-1i*H*tstep);

    ctemp = c0.';
    cs = zeros(length(ts),length(c0));

    cs(1,:) = c0.';

    for j=1:length(ts)-1

        ctemp = Ustep*ctemp;
        cs(j+1,:) = ctemp.';
    end

end