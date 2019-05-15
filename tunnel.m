function tunnel()

    % performance
    dx = 0.05;
    xmax = 10;
    dt = 1;   
    nstep = 300;
    dstep = 1;

    % particle info
    s = 1;
    p0 = 20;
    m = 100;
    
    % Barrier is around X=0
    L = 10;
    V0 = 0.6;

    [x, psi] = gaussianwave(xmax, dx, s, p0, -5);
    
    N = length(x);
    xleft = x(1:round(N/2));
    xright = x(round(N/2)+1:N);

    % We already know the eigenstates:

    p = linspace(0, 2*pi/dx, N) - (pi/dx - 2*pi/dx/N);
    dp = p(2) - p(1);
    p1 = sqrt(p.^2 - 2*m*V0);
    r = (p.^2-p1.^2).*(1-exp(2i*p1*L))./((p+p1).^2 - (p-p1).^2.*exp(2i*p1*L));
    t = 4*p.*p1.*exp(-1i*(p-p1)*L)./((p+p1).^2 - (p-p1).^2.*exp(2i*p1*L));

    % row: different p; column: different x
    basis = [exp(1i*(p'*xleft)) + diag(r)*exp(-1i*(p'*xleft)), diag(t)*exp(1i*(p'*xright))];
    
    % normalization
    for j = 1:N
        basis(j,:) = basis(j,:) / norm(basis(j,:)) / sqrt(dx);
    end
    
    % a0: Using inner product
    % a0 = <p|\psi>
    a0 = conj(basis)*transpose(psi)*dx;

    % Check if it is right
    % |psi> = a0|p>
    hold on;
    plot(x, abs(psi), 'r-')
    plot(x, abs(transpose(a0)*basis), 'b-');
    hold off;

    % Evolve
    for j = 1:dstep:nstep
        E_coeff = diag(exp(-1i*p.^2/2/m*dt*j));
        plot([0,0], [0,1], 'k-', x, abs(transpose(a0)*(E_coeff*basis)));
        ylim([0,1])
        pause(0.1);
    end

end