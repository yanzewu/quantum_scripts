function [ts, as] = runao(H, a0, t, tstep, tstepc)
% Run a system of harmonic oscillators by Schodinger's equation
% da/dt = -i*H*a;
% Args:
% a0: NxN matrix, a(i,j) = <a0*(i)a(j)>
% t: is total time;
% tstep: Output step;
% tstepc (optional): Step of ODE solving, must be integer multiples of tstep.
%
% Returns are similar to ode() in matlab;
% ts: List of time;
% cs: List of a at different times, with each row as a(:)

    if ~exist('tstepc', 'var')
        tstepc = t;
    end

    as = [];
    ts = [];
    a0temp = a0(:);
    t0temp = 0;

    sz = size(a0);

    while t0temp < t
        [ttemp, atemp] = run_session(H, reshape(a0temp, sz), min(t-t0temp, tstepc), tstep);
        as = [as;atemp(1:end-1,:)];
        ts = [ts;ttemp(1:end-1)+t0temp];
        a0temp = atemp(end,:);
        t0temp = ttemp(end)+t0temp;
    end

    as = [as; a0temp];
    ts = [ts; t0temp];
    
end

function [t_, a] = run_session(H, a0, t, tstep)

    sz = size(a0);
    col = sz(1);

    odeopt = odeset('RelTol', 1e-8, 'AbsTol', 1e-9);

        function dy = evolve(~, y, H_, col_)
            y_ = reshape(y, [col_,col_]);
            dy = zeros(col_,col_);
            for j = 1:col_
                dy(j,:) = -1i*y_(j,:)*H_;
            end
            dy = dy(:);
        end

    [t_, a] = ode113(@(~,y) evolve(0,y,H,col), 0:tstep:t, a0(:), odeopt);

end