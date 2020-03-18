% construct plant-taper window
% from 0 to N-1
function w = planktaper(N,eps)
w = zeros(1,N);
w(1) = 0;
for n=1:(floor(eps*(N-1))-1)
	Zp = 2*eps*( 1/(1+ (2*n/(N-1)-1)) + 1/(1 - 2*eps+ ( 2*n/(N-1)-1)) );
	w(n+1) = 1/(exp(Zp)+1);
end

for n=(floor(eps*(N-1))):(floor((1-eps)*(N-1)))
	w(n+1) = 1;
end
for n = (floor((1-eps)*(N-1))+1):(N-2)
	Zm = 2*eps*( 1/(1- (2*n/(N-1)-1)) + 1/( 1 - 2*eps- ( 2*n/(N-1)-1)) );
	w(n+1) = 1/(exp(Zm)+1);
end
w(N) = 0;
