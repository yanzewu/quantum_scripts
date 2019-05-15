function E = get_env(d, N)
% Set up environment around E=0;
% d: Energy interval;
% N: Number of state, must be odd;

N1 = round((N-1)/2);
E = -N1*d:d:N1*d;

end