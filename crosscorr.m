% circular autocorrelation
% need f and g to be of the same length
function z = crosscorr(f,g)
if (length(f) ~= length(g))
	disp('vector lengths do not match');
else
	N = length(f);
	z = zeros(1,N);
	for i=1:N
	z(i) = dot(f,shift(g,i));
	end
end
