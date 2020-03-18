% compute the cumulative power distribution
function s = cupow(z) 
s = zeros(1,length(z));
for i=1:length(z)
	s(i) = sum(abs(z(1:i)));
end


