function wavepacketms()
% Multi-surface evolution

    % performance
    dx = 0.02;
    xmax = 15;
    dt = 0.01;
    nstep = 5000;    % total steps
    dstep = 10;     % display steps
    
    % particle info
    s = 1;
    p0 = 50;
    m = 100;

    [x, psi] = gaussianwave(xmax, dx, s, p0, -5);
    psi = psi./(norm(psi)*sqrt(dx));
    
    N = length(x);
   
    % Diabatic Hamiltonian
    K = -get_laplacian(1, 2)./(2*m*dx^2);
    V = [(x+5).^2*0.02;(x-5).^2*0.02]*0.6;  % diabatic surface
    S = [0,0.1;0.1,0];  % diabatic coupling

    nsurf = length(S);
  
    % diabatic evolution
    function dY = evolve_d(~, Yarr, K_, V_, S_, sz_)
        Ymat = reshape(Yarr, sz_);
        
        dY = -1i*(conv2(Ymat, K_, 'same') + V_.*Ymat + S_*Ymat);
        dY = dY(:);
    end

    
    % Adiabatic surface and eigenvectors
    V_ad = zeros(nsurf, N);
    basis_ad0 = zeros(nsurf, N);
    basis_ad1 = zeros(nsurf, N);
    basis_tad0 = zeros(nsurf, N);
    basis_tad1 = zeros(nsurf, N);
    
    for i = 1:1:N
        V_local = S + diag(V(:,i));
        [evec, evalue] = sorted_eig(V_local);
        V_ad(:,i) = evalue;
        basis_ad0(:,i) = evec(:,1);
        basis_ad1(:,i) = evec(:,2);
        basis_tad0(:,i) = evec(1,:)';
        basis_tad1(:,i) = evec(2,:)';
    end

    % Derivative Coupling
    dR = get_fs_mat(0,1,4)/dx;
    D = [dot(basis_ad1,conv2(basis_ad0, dR, 'same')); dot(basis_ad0,conv2(basis_ad1, dR, 'same'))];
    
    % adiabatic evolution
    function dY = evolve_ad(~, Yarr, K_, V_, D_, DR_, sz_)
        Ymat = reshape(Yarr, sz_);
        
        DC = D_.*conv2(Ymat, DR_, 'same');
        dY = -1i*(conv2(Ymat, K_, 'same') + V_.*Ymat + DC([2,1],:));
        dY = dY(:);
    end

    
    % Initialize state   
    state = zeros(nsurf, 1);
    state(1) = 1;   % initial state
    psi_d = state*psi;
    
    % Re-projection
    psi_ad = basis_tad0.*psi_d(1,:) + basis_tad1.*psi_d(2,:);

    sz_psi = size(psi_d);
    evolve_d_dummy = @(t, Yarr) evolve_d(t, Yarr, K, V, S, sz_psi);
    evolve_ad_dummy = @(t, Yarr) evolve_ad(t, Yarr, K, V, D, dR, sz_psi);
    
    for i = 1:1:nstep
        
        [~, psiarray] = ode23(evolve_d_dummy, [0, dt], psi_d(:));
        psi_d = reshape(psiarray(size(psiarray,1),:), sz_psi);
        
        [~, psiarray] = ode23(evolve_ad_dummy, [0, dt], psi_ad(:));
        psi_ad = reshape(psiarray(size(psiarray,1),:), sz_psi);
        
        if mod(i, dstep) == 0
            hold off
            
            % diabatic states
            plot(x, abs(psi_d(1,:)), 'r-', x, abs(psi_d(2,:)), 'b-');
            hold on
            
            % potential surfaces
            plot(x, V(1,:), 'k-', x, V(2,:), 'k-');
            plot(x, V_ad(1, :), 'k--', x, V_ad(2,:), 'k--');

            % adiabatic states projection
            psi_d_recons = basis_ad0.*psi_ad(1,:) + basis_ad1.*psi_ad(2,:);
            plot(x, abs(psi_d_recons(1,:)), 'r--', x, abs(psi_d_recons(2,:)), 'b--');
            
            ylim([0,1]);
            
            % Additional Normalization
            psi_d = psi_d ./ (normvec(psi_d)*sqrt(dx));
            pause(0.1);
        end
    end
    
end

function nv = normvec(v)
% sqrt(sum(abs(v).^2))

    v2 = abs(v).^2;
    nv = sqrt(sum(v2(:)));
end

function [V,D] = sorted_eig(m)

    [v, d] = eig(m);
    [D, p] = sort(diag(d));
    V = v(:,p);
end