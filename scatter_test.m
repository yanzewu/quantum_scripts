
function scatter_test()

    L = 0.1;
    V0 = 1;

    V = @(x) and(x>=0, x<L).*V0;
    
    N = 200;
    k = 2;
    E = k*k/2;
    [refl, trans, psi] = m_scatter(N, k, V);
    trans
    m_trans = 1/(1+V0^2*sin(sqrt(2*E-V0)*L)^2/4/E/(E-V0))
    trans + refl
    plot(real(psi))

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
        H(1,q+1) = sum(exp(-1i*k*dx*((-n+1:0)-q)).*get_t(q-(-n+1:0)));
    end
    for q = N+1-max_stencil:N+1
        H(N+2,q+1) = sum(exp(1i*k*dx*((0:n-1)+q)).*get_t(N+1-q+(0:n-1)));
    end
    
    % LHS: Inside stencils
    for q = 2:N+1
        reallimit = min(q+max_stencil-1,N+1);
        H(q,q:reallimit) = get_t((q:reallimit)-q);
    end
    
    H = H + H' - diag(diag(H)');
    H(2:N+1,2:N+1) = H(2:N+1,2:N+1) + diag(V(((2:N+1)-(N+3)/2)*dx));
    
    % LHS: Energy
    lhs = H - k*k/2/m*eye(N+2);
    
    % RHS
    rhs = zeros(N+2,1);
    for q = 0:max_stencil-1
        rhs(q+1) = sum(exp(1i*k*dx*((-n+1:0)-q)).*get_t(q-(-n+1:0)));
    end
    rhs(1) = rhs(1) + k*k/2/m;
    
    % Calculation
    ret = lhs \ rhs;
    refl = abs(ret(1)).^2;
    trans = abs(ret(N+2)).^2;
    psi = ret(2:N+1);

end
