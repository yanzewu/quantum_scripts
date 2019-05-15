function [ts, cs] = runho(H, c0, t, tstep, tstepc)
% Run a Hamiltonian by Schrodinger's Equation:
% idc/dt = Hc;
% Args:
% c0: Initial vector of c;
% t: Total time;
% tstep: Output step;
% tstepc (optional): Step of ODE solving, must be integer multiples of tstep.
%
% Returns are similar to ode() in matlab;
% ts: List of time;
% cs: List of c at different times, with each row as c(t);

    if ~exist('tstepc', 'var')
        tstepc = t;
    end

    cs = [];
    ts = [];
    c0temp = c0;
    t0temp = 0;

    while t0temp < t
        [ttemp, ctemp] = run_session(H, c0temp, min(t-t0temp, tstepc), tstep);
        cs = [cs;ctemp(1:end-1,:)];
        ts = [ts;ttemp(1:end-1)+t0temp];
        c0temp = ctemp(end,:);
        t0temp = ttemp(end)+t0temp;
    end

    cs = [cs; c0temp];
    ts = [ts; t0temp];

end

function [t_, c] = run_session(H, c0, t, tstep)

    odeopt = odeset('RelTol', 1e-9, 'AbsTol', 1e-10);
    [t_, c] = ode45(@(~,y) -1i*H*y, 0:tstep:t, c0, odeopt);

end