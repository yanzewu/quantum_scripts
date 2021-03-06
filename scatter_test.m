
function scatter_test()

    L = 0.1;
    V0 = 1;

    V = @(x) and(x>=0, x<L).*V0;
    
    N = 200;
    
    kall = 0.1:0.2:2;
    transall = [];
    m_transall = [];
    for k = kall
        E = k*k/2;
        [refl, trans, psi] = m_scatter(N, k, V);
        m_trans = 1/(1+V0^2*sin(sqrt(2*(E-V0))*L)^2/4/E/(E-V0));
        transall = [transall;trans];
        m_transall = [m_transall; m_trans];
    end
    plot(kall, transall, 'rx', kall, m_transall, 'b-');

end

function [refl, trans, psi] = m_scatter(N, k, V)
% Solve scatter equation
% (H-E)\psi = (E-H)\psi_{inc}

    n = 3;
    m = 1;
    dx = 0.02;
    
    t = -1/2/m/dx/dx*[-5/2,4/3,-1/12];  % t(n+1) -> T(a,a+n)=T(a,a-n)
    
    H = zeros(N+2);
    
    function m_t = get_t(d)
        m_t = zeros(1, length(d));
        for j = 1:length(d)
            if abs(d(j)) <= 2
                m_t(j) = t(abs(d(j))+1);
            end
        end
    end
    
    max_stencil = 3;

    % LHS: Outside stencils
    for q = 0:max_stencil-1
        H(q+1,1) = sum(exp(-1i*k*dx*(-n+1:0)).*get_t(q-(-n+1:0)));
    end
    for q = N+1-max_stencil:N+1
        H(q+1,N+2) = sum(exp(1i*k*dx*(0:n-1)).*get_t(N+1-q+(0:n-1)));
    end
    
    % LHS: Inside stencils
    for q = 1:N+2
        min_q = max(q-max_stencil+1, 2);
        max_q = min(q+max_stencil-1, N+1);
        H(q,min_q:max_q) = get_t((min_q:max_q)-q);
    end
    
    H(2:N+1,2:N+1) = H(2:N+1,2:N+1) + diag(V(((2:N+1)-(N+3)/2)*dx));
    
    % LHS: Energy
    lhs = H - k*k/2/m*eye(N+2);
    
    % RHS
    rhs = zeros(N+2,1);
    for q = 0:max_stencil-1
        rhs(q+1) = -sum(exp(1i*k*dx*(-n+1:0)).*get_t(q-(-n+1:0)));
    end
    rhs(1) = rhs(1) + k*k/2/m;
    
    % Calculation
    ret = lhs \ rhs;
    refl = abs(ret(1)).^2;
    trans = abs(ret(N+2)).^2;
    psi = ret(2:N+1);

end
