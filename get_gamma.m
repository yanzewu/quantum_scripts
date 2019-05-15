function [Gamma, dE] = get_gamma(E0, V, E)
% Calculate gamma and energy shift of relaxation by
% Gamma = 2*pi*|Vj|^2, Ej=E;
% dE = sum_{k/=j}{|Vk|^2/(Ek-E)}, Ej=E;
%
% Args:
% E0: Energy of single state;
% V: Coupling constant, can be either a number or a list with size same as E;
% E: List of energy, must be monotonic increasing;

    E0_idx = find(E>E0, 1)-1;
    d0 = E(E0_idx+1) - E(E0_idx);
    
    if length(V) == 1
        Gamma = 2*pi*abs(V)^2/d0;
        if E(E0_idx) == E0
            dE = (sum(1./(E0-E(1:E0_idx-1)))+sum(1./(E0-E(E0_idx+1:end))))*abs(V)^2;
        else
            dE = sum(1./(E0-E))*abs(V)^2;
        end
    else
        Gamma = 2*pi*abs(V(E0_idx))^2/d0;
        if E(E0_idx) == E0
            dE = (sum(abs(V).^2./(E0-E(1:E0_idx-1)))+sum(abs(V).^2./(E0-E(E0_idx+1:end))));
        else
            dE = sum(abs(V).^2./(E0-E));
        end
    end
    
    
end