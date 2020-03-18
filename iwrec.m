% recursive implementation of inverse wavelet transform
% wavelet transform of z wrt parent wavelets u and v with smallest possible dimension sdim
% e.g. for p-th stage wavelets, sdim = N/2^(p-1)
function w = iwrec(z, sdim, u, v)
n = length(z);
% first part of the sum
first = realconv(up(z(1:n/2)),v);
if n <= 2*sdim
	% stopping condition
	second = realconv(up(z(n/2+1:n)),u);
else
	% recursion
	second =realconv(up(iwrec(z(n/2+1:n), sdim, fold(u), fold(v))), u);
end
w = first + second;
