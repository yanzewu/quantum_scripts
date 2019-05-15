function model1a()

    V = 0.1;
    d = 0.1;
    N = 201;
    randV = true;
    colors = ones(length(V),1)*'r';

    hold off
    cla();
    box on;
    
    for j=1:length(V)
        V_ = V(j);
        d_ = d(j);
        [t,p1,penv,gamma] = calc_model(20, 0.1, V_, d_, N, randV);
        gamma
        figure(1)
        hold on
        plot(t,p1, strcat(colors(j),'-'), 'LineWidth',2);
        plot(t,exp(-gamma*t),'k--','LineWidth',1)
    end
    xlabel('t')
    ylabel('Probability')
    set(gca, 'FontSize', 14)

end

function [t,p1,penv,gamma] = calc_model(time, tstep, V, d, N, randV)

    E = get_env(d,N);
    [gamma,~] = get_gamma(0,V,E);
    E = E + rand(1,N)*0.1;
    H = zeros(N+1,N+1);
    H(1,1) = 0;
    
    if randV
        H(1,2:end) = (2*round(rand(1,N))-1)*V;
        H(2:end,1) = H(1,2:end);
    else
        H(1,2:end) = V;
        H(2:end,1) = V; 
    end

    H(2:end,2:end) = diag(E);
    [t,c] = runhm(H, [1,zeros(1,N)], time, tstep);
    p1 = abs(c(:,1)).^2;
    penv = sum(abs(c(:,2:end)).^2, 2);
    return

end