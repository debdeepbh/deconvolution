% crude derivative
function w = deriv(z)
N = length(z);
k=1;
for i=1:N-k
	w(i) = (z(i+k) - z(i))/k;
end

