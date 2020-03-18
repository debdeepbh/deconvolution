% given the mother(or father) wavelet get the other
function v = getother(u) 
N = length(u);
% v(k) = (-1)^{k-1} u_bar(1-k)
for k=0:N-1
	% mod keep the value in the range of 0 and N-1
	% +1 shifts the vector so that u is defined
	% place the result in v(k+1)
	v(k+1) = (-1)^(k-1)*conj(u(mod(1-k,N) +1));
end

