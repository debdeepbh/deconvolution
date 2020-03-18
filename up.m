% upsampling of a vector of even dimension
function Uz = up(z)
n = length(z);
for k = 1:n
	Uz(2*k -1) = z(k);
	Uz(2*k) = 0;
end
