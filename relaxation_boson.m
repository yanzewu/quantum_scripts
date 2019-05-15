function relaxation_boson

    w0 = 20;
    N0 = 5;
    V = 0.1;
    d = 0.1;
    T = 50;
    N = 101;
    we_shift = 0.0;

    [w0e, n0e] = get_env_boson(w0+we_shift, d, N, 1/T);

    [gamma, dw0] = get_gamma(w0, V, w0e)

    H = zeros(N+1);
    H(1,1) = w0;
    H(1,2:end) = V;
    H(2:end,1) = V;
    H(2:end,2:end) = diag(w0e);

    [ts, as] = runam(H, eye(N+1), 60, 1);

    nvec = [N0,n0e];

    a1s = as(:,1:N+1);
    n1s = (abs(a1s).^2*(nvec)');

    n1t = N0*exp(-gamma*ts)+(1-exp(-gamma*ts))/(exp((w0+dw0)/T)-1);
    n1t2 = exp(-gamma*ts)*N0;
    for j = 1:N
        n1t2 = n1t2 + V^2/((w0+dw0-w0e(j))^2+gamma^2/4)*n0e(j).*(1+exp(-gamma*ts)-2*exp(-0.5*gamma*ts).*cos((w0+dw0-w0e(j))*ts));
    end

    figure(1)
    plot(ts, n1t, 'k--', ts, n1t2, 'b--', ts, n1s, 'r-');
    xlabel('t')
    ylabel('Population')
    legend('Prediction 1', 'Prediction 2', 'Numerical')

    effw = T*log(1./n1s+1);
    figure(2);
    plot(ts, ones(size(ts))*(w0+dw0), 'k--', ts, effw, 'r-');
    xlabel('t')
    ylabel('Effective energy')

end