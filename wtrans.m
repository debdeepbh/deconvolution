%wavelet transform, general
% Caution: z should be a row vector, not column vector.
% size of N should be divisible by 2^(p-1)
function w = wtrans(z, type, p)
N = length(z);
[u, v] = filt(type, N);

util = ifft(conj(fft(u)));
vtil = ifft(conj(fft(v)));


% computing the smallest dimension based on p
w = wrec(z, N/2^(p), util, vtil);

