% downsampling of a vector of even dimension
function Dz = down(z)
% even check
if mod(length(z),2) == 1
	printf 'error:odd length\n'
else
	n = length(z)/2;
	for k = 1:n
		Dz(k) = z(2*k - 1);
	end
end
