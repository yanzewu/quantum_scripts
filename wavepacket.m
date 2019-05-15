function [T,KE,PE,Etot] = wavepacket()

    dx = 0.05;
    xmax = 10;
    dt = 1;
    nstep = 3000;    % total steps
    dstep = 10;     % display steps
    cstep = 500;   % computation steps
    acc = 4;
    
    % particle info
    s = 1;
    p0 = 20;
    m = 100;

    [x, psi] = gaussianwave(xmax, dx, s, p0, -8);
    psi = psi./(norm(psi)*sqrt(dx));

    N = length(x);
   
    % Create Hamiltonian
    K = -get_fs_mat(N, 2, acc, true)./(2*m*dx^2);
    K_kernel = -get_fs_mat(N, 2, acc, false)./(2*m*dx^2);
    %V0 = get_potential(x, 'harmonic', [0.25, 0])';
    V0 = get_potential(x, 'none', []);

    V0(1:0.25*N) = -1.9;
    V0(0.75*N:N) = -2.1;
    
    V0 = V0.';
    V = diag(V0);
    
    % Solving ODE
    function yp = evolve(~, y)
       yp = -1i*(conv(y, K_kernel, 'same') + V0.*y);
    end

    KE = zeros(1, nstep+1);
    PE = zeros(1, nstep+1);
    Etot = zeros(1, nstep+1);
    T = linspace(0, nstep, nstep+1) * dt;

    odeopt = odeset('RelTol', 1e-5);
    
    psi = transpose(psi);
    
    for n = 0:cstep:nstep-1
      
        lenret = min(nstep-n, cstep)+1;
        [~, ret] = ode23(@evolve, linspace(0, lenret-1, lenret)*dt, psi, odeopt);
        
        psi = transpose(ret(lenret, :));
        psi = psi / norm(psi) / sqrt(dx);
        
        % Energy Calculation
        for j = 1:lenret
            KE(n + j) = real(ret(j,:)*K*ret(j,:)'*dx);
            PE(n + j) = real(ret(j,:)*V*ret(j,:)'*dx);
        end
        
        Etot = KE + PE;

        % Animation
        for j = 0:dstep:lenret-1
            fprintf('t=%.4f\n', (n+j)*dt);

            subplot(2,1,1);
            cla();
            hold on
            plot(T(1:n+j+1),KE(1:n+j+1), 'r-');
            plot(T(1:n+j+1),PE(1:n+j+1), 'b-');
            plot(T(1:n+j+1),Etot(1:n+j+1) ,'m-');
            legend('KE','PE','Etot');
            hold off

            subplot(2,1,2);
            cla();
            hold on
            plot(x, (V0-min(V0))/(max(V0)-min(V0)), 'r-');
            plot(x, abs(ret(j+1,:)), 'k-');
            ylim([0 2]);
            hold off

            drawnow;
            pause(0.1);
        end
    
    end
end