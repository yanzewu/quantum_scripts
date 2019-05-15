function [ts, as] = runam(H, a0, t, tstep)
% Run a system of harmonic oscillators by Schodinger's equation
% da/dt = -i*H*a;
% Args:
% a0: NxN matrix, a(i,j) = <a0*(i)a(j)>
% t: is total time;
% tstep: Output step;
%
% Returns are similar to ode() in matlab;
% ts: List of time;
% cs: List of a at different times, with each row as a(:)

    if ~exist('tstep', 'var')
        tstep = t;
    end

    ts = 0:tstep:t;
    as = zeros(length(ts),numel(a0));

    Ustep = expm(-1i*H*tstep).';
    atemp = a0;

    as(1,:) = a0(:);

    for j=1:length(ts)-1
        atemp = eval_a(Ustep, atemp);
        as(j+1,:) = atemp(:);
    end

end

function a = eval_a(Ut, a0)

    a = zeros(size(a0));
    for j=1:length(a0)
        a(j,:) = a0(j,:)*Ut;
    end

end