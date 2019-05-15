
% create a function exp(-x^2/s^2+ikx), normed.
function [x, psi] = gaussianwave(xmax, dx, sigma, k, center)

    if ~exist('center', 'var')
        center = 0;
    end
    
    % center
    N = round(xmax * 2 / dx) + 1;
    x = linspace(-xmax, xmax, N);
    psi = exp(-(x-center).^2/2/sigma^2 + 1i*k*(x-center))/(sigma^0.5*pi^0.25);

end