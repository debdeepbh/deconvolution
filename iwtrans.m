% inverse wavelet transform, general
function w = iwtrans(z, type, p)
N = length(z);
[u, v] = filt(type, N);
w = iwrec(z, N/2^(p), u, v);

