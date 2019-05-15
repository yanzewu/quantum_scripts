function [w, n0] = get_env_boson(w0, d, N, beta)

N1 = round((N-1)/2);
w = w0-N1*d:d:w0+N1*d;
n0 = 1./(exp(beta*w)-1);

end