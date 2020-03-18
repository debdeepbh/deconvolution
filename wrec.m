% recursive implementation of wavelet transform
% wavelet transform of z wrt parent wavelets u and v with smallest possible dimension sdim
% e.g. for p-th stage wavelets, sdim = N/2^(p-1)
function w = wrec(z, sdim, util, vtil)
first = down(realconv(z,vtil));

if length(z) <= 2*sdim
	% some stopping condition
       	second = down(realconv(z, util));
	length(second);
else
	% recursion
       	second = wrec(down(realconv(z,util)), sdim, fold(util),fold(vtil));
end




w = [first second];

