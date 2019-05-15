function wavepacket2d()

    % Performance

    xmax = 2;
    dx = 0.05;
    dt = 1;
    nstep = 1000;
    dstep = 1;
    
    % Motion
       
    v = [0,0];
    center = [0,0];
    m = 40;
    mz = 2;
    k = 4;
    sigma = [m*k,m*k].^(-0.25);
    
    % Code
    
    [x, y, psi] = gaussianwave2d(xmax, xmax, dx, sigma, v.*m, center);
    angle = atan2(y, x);
    r = sqrt(x.^2+y.^2);
    
    psi = psi .* exp(1i*mz*angle) .* r.^mz;
 %   psi = psi .* ((2*y^2-1).*(2*x.^2-1) + 50*exp(1i*6*angle) .* r.^6);
    psi = psi ./ (normvec(psi)*dx);
    
    
    sz_psi = size(psi);
    
    KE0 = norm(v).^2*m/2
    
    K = -get_laplacian(2, 2)/dx^2/2/m;
    V = 0.5 * k * (x.^2+y.^2);
    
    function dY = evolve(t, Yarr, K_, V_, sz_)
        Ymat = reshape(Yarr, sz_);
        dY = -1i*(conv2(Ymat, K_, 'same') + Ymat.*V_);
        dY = dY(:);
    end
    
    evolve_dummy = @(t, Yarr) evolve(t, Yarr, K, V, sz_psi);

    mesh(x, y, abs(psi));
    
    for i = 1:1:nstep
        
        psiflat = psi(:);
        
        [~, psiarray] = ode23(evolve_dummy, [0, dt], psiflat);
        [rows_, ~] = size(psiarray);
        psiflat = psiarray(rows_,:);
        
        psi = reshape(psiflat, sz_psi);
        
        if mod(i, dstep) == 0
            mesh(x, y, real(psi));
            zlim([-2,2]);
 %           imagesc(abs(psi));
            psi = psi ./ (normvec(psi)*dx);
            pause(0.1);
        end
    end

end

function [X1, X2, psi] = gaussianwave2d(xmax1, xmax2, dx, sigma, k, center)
% sigma, k, center is required to be 1x2
    
    N1 = round(xmax1 * 2 / dx) + 1;
    N2 = round(xmax2 * 2 / dx) + 1;
    x1 = linspace(-xmax1, xmax1, N1);
    x2 = linspace(-xmax2, xmax2, N2);
    
    [X1, X2] = meshgrid(x1, x2);

    psi = exp(-(X1-center(1)).^2/2/sigma(1)^2 + 1i*k(1)*(X1-center(1))-(X2-center(2)).^2/2/sigma(2)^2 + 1i*k(2)*(X2-center(2)));
    
    psi = psi ./ (normvec(psi)*dx);
    
end

function nv = normvec(v)

    v2 = abs(v).^2;
    nv = sqrt(sum(v2(:)));
end