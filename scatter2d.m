function scatter2d

    L = 0.1;
    V0 = 1;
    
    Vx = @(x) and(x>=0, x<L).*V0;
    Vy = @(y) zeros(size(y));
    
    Nx = 100;
    Ny = 20;
    Ly = 0.02*(Ny-1);
    
    k = 2;
    
    kall = 0.1:0.5:2;
    transall = [];
    m_transall = [];
    for k = kall
    
    E = k*k/2;
    
    Ey = pi^2/2/Ly^2;
    y = (0:Ny-1)*0.02;
    psi0_y = sin(pi/Ly*y) * sqrt(2/Ly);
    
    [refl, trans, psi] = m_scatter(Nx, Ny, k, Ey, psi0_y.', Vx, Vy);
    trans
    m_trans = 1/(1+V0^2*sin(sqrt(2*(E-V0))*L)^2/4/E/(E-V0))
    trans+refl
    transall = [transall; trans];
    m_transall = [m_transall; m_trans];
    end
    plot(kall, transall, 'rx', kall, m_transall, 'b-');
    %surf(abs(psi).^2);

end

function [refl, trans, psi] = m_scatter(Nx, Ny, k, Ey, psi0_y, Vx, Vy)
% Pseudo 2d scattering
% psi0_y is required as column vector

    n = 3;
    m = 1;
    dx = 0.02;
    x = dx*((1:Nx) - (Nx+1)/2);
    y = dx*(0:Ny-1);
    
    % bound states
    Hy = -1/2/m/dx/dx*get_fs_mat(Ny, 2, 4, true);
    [V__, D__ ] = eigs(Hy, 1, 'smallestreal');
    Ey = D__(1,1);
    psi0_y = V__(:,1);
    %figure(2);
    %plot(real(psi0_y));
    
    t = -1/2/m/dx/dx*[-5/2, 4/3, -1/12]; % 1d
    function m_t = get_t(d)
        m_t = zeros(size(d));
        for j = 1:length(d)
            if abs(d(j)) <= 1
                m_t(j) = t(abs(d(j))+1);
            end
        end
    end

    t2 = zeros(length(t));
    t2(:,1) = t;
    t2(1,:) = t2(1,:) + t;
    function m_t2 = get_t2(del_x, del_y)
        % return matrix with shape (length(dx),length(dy))
        m_t2 = zeros(length(del_x), length(del_y));
        for j = 1:length(del_x)
            for j2 = 1:length(del_y)
                if and(abs(del_x(j)) <= 2, abs(del_y(j2)) <= 2)
                    m_t2(j, j2) = t2(abs(del_x(j))+1, abs(del_y(j2)) + 1);
                end
            end
        end
    end
    
    N = Nx * Ny;
    H = zeros(Nx*Ny + 2);
    max_stencil = 3;
    
    gr_x = -n+1:0;
    gt_x = 0:n-1;
    psii_x = exp(1i*k*dx*gr_x);
    psir_x = exp(-1i*k*dx*gr_x);
    psit_x = exp(1i*k*dx*gt_x);
    
    gr_y = 1:Ny;
    gt_y = 1:Ny;
    psii_y = psi0_y;
    psir_y = psi0_y;
    psit_y = psi0_y;
    
    % Corners 
    % H(1,1) = Ey + sum(psir_x .* get_t(gr_x));
    % H(N+2,N+2) = Ey + sum(psit_x .* get_t(gt_x));
    H(1,1) = sum(Vy(y).'.*abs(psir_y).^2);
    for qy = 1:Ny
        for qy2 = 1:Ny
            H(1,1) = H(1,1) + sum(psir_x.* get_t2(gr_x, qy-qy2).') * conj(psir_y(qy))*psir_y(qy2);
        end
    end
    H(N+2,N+2) = sum(Vy(y).'.*abs(psit_y).^2);
    for qy = 1:Ny
        for qy2 = 1:Ny
            H(N+2,N+2) = H(N+2,N+2) + sum(psit_x.*get_t2(gt_x, qy-qy2).') * conj(psit_y(qy))*psit_y(qy2);
        end
    end
    
    loc2 = @(a, b) (a-1)*Ny + b;    % maps 2D coords to 1D
    
    % Edges
    for qx = 1:max_stencil-1   % <x|H|\psi_r>
        for qy = 1:Ny
            H(loc2(qx,qy) + 1, 1) = psir_x * get_t2(qx-gr_x, qy-gr_y) * psir_y;
        end
    end
    for qx = Nx-max_stencil+1:Nx  % <x|H|\psi_t>
        for qy = 1:Ny
            H(loc2(qx,qy) + 1, N+2) = psit_x * get_t2(Nx+1-qx+gt_x, qy-gt_y) * psit_y;
        end
    end
    for qx = 1:max_stencil-1  % <xr|H|x>
        for qy = 1:Ny
            H(1, loc2(qx,qy) + 1) = sum(get_t2(qx, qy-gr_y)*conj(psir_y));
        end
    end
    for qx = Nx-max_stencil+1:Nx % <xt|H|x>
        for qy = 1:Ny
            H(N+2, loc2(qx,qy) + 1) = sum(get_t2(Nx+1-qx, qy-gt_y)*conj(psit_y));
        end
    end
    
    % Inside
    for qx = 1:Nx
        for qx2 = max(1, qx-max_stencil+1):min(Nx, qx+max_stencil-1)
            for qy = 1:Ny
                min_qy2 = max(1, qy-max_stencil+1);
                max_qy2 = min(Ny, qy+max_stencil-1);
                H(loc2(qx, qy) + 1, loc2(qx2, min_qy2)+1:loc2(qx2, max_qy2)+1) = get_t2(qx2-qx, qy-(min_qy2:max_qy2));
            end
        end
    end
    
    % Potential

    H(2:N+1, 2:N+1) = H(2:N+1, 2:N+1) + kron(diag(Vx(x)), eye(Ny)) + kron(eye(Nx), diag(Vy(y)));
    
    % Energy
    lhs = H - (k*k/2/m+Ey)*eye(N+2);
    
    % RHS
    rhs = zeros(N+2, 1);
    %rhs(1) = -(Ey + sum(psii_x .* get_t(gr_x)));
    rhs(1) = -sum(Vy(y).'.*abs(psii_y).^2);
    for qy = 1:Ny
        for qy2 = 1:Ny
            rhs(1) = rhs(1) - sum(psii_x.*get_t2(gr_x, qy-qy2).') * conj(psir_y(qy))*psii_y(qy2);
        end
    end
    for qx = 1:max_stencil - 1
        for qy = 1:Ny
            rhs(loc2(qx,qy) + 1) = -psii_x * get_t2(qx-gr_x, qy-gr_y) * psii_y;
        end
    end
    rhs(1) = rhs(1) + k*k/2/m + Ey;
    
    % Calculation
    ret = lhs \ rhs;
    refl = abs(ret(1)).^2;
    trans = abs(ret(N+2)).^2;
    psi = reshape(ret(2:N+1), [Ny, Nx]);

end