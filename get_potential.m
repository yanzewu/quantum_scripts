function v = get_potential(x, type, args)
% x: Array of x coordinates
% type: Potential type, 'none'/'harmonic'/'morse'
% args: none==>[]; harmonic==>[k,c]; morse==>[d, alpha]; barrier==>[V0,xl,xr]


    switch type
        case 'none'
            v = zeros(1, length(x));
        case 'harmonic'
            v = 0.5*args(1)*(x-args(2)).^2;
        case 'morse'
            v = args(1)*(1-exp(-args(2)*x)).^2;
        case 'barrier'
            v = zeros(1, length(x));
            left = round((args(2) - x(1))/(x(2) - x(1))) + 1;
            right = round((args(3) - x(1))/(x(2) - x(1))) + 1;
            v(left:right) = args(1);
        otherwise
            error('Undefined potential');
    end

end