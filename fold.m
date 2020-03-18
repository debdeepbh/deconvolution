% folding a vector at half
function y = fold(z)
n = length(z);
% Check if even
if mod(n,2) != 0
	printf 'error: vector not even';
else
	a = z(1:n/2);
	b = z(n/2+1:n);
	y = a + b;
end



