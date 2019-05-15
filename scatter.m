function [refl, trans] = scatter(E)
    % Computation of scattering state

    % Performance
    dx = 0.05;
    xmax = 10;
    acc = 4;
    ext = 5; % extra grid number in each side
    ext_est = 2;    % extra grid number in estimation. should < ext/2

    % particle info
    m = 10;

    x = -xmax-ext*dx:dx:xmax+ext*dx;
    N = length(x);
    N_reg = length(x)-ext*2;

    % Extra matrix
    k = sqrt(E*2*m);
    basis_inc = [exp(1i*k*x(1:ext)),zeros(1,N-ext)];
    basis_r = exp(-1i*k*x(1:ext));
    basis_t = exp(1i*k*x(N-ext+1:N));

    M_ext = [zeros(ext,N_reg),basis_r.',zeros(ext,1); eye(N_reg),zeros(N_reg,2); zeros(ext,N_reg+1),basis_t.'];

    % H
    K = -get_fs_mat(N, 2, acc, true)./(2*m*dx^2);
    V0 = get_potential(x, 'barrier', [2, 0, 1]);
    %V0 = exp(-(x-2).^2)*2+exp(-(x+2).^2)*2+exp(-(x-6).^2)*2;% An interesting barrier
    H = K + diag(V0);

    % H real
    H_left = (H - eye(N)*E)*M_ext;
    right = -(H - eye(N)*E)*basis_inc.';

    % 
    ret = H_left(ext-ext_est:N-ext+1+ext_est,:) \ right(ext-ext_est:N-ext+1+ext_est);
 
    refl = abs(ret(N_reg+1))^2;
    trans = abs(ret(N_reg+2))^2;
    
    psi0 = M_ext*ret + basis_inc.';
    
    for t = 0:0.1:10
        psi = psi0 * exp(-1i*E*t);
        plot(x, (V0-V0(1))/(max(V0)-min(V0)), 'k-');
        hold on;
        plot(x, real(psi), 'b-');
        plot(x, imag(psi), 'r-');
        ylim([-2,2]);
        legend('V', 'Re(\psi)', 'Im(\psi)')
        hold off;
        pause(0.1);
    end

end