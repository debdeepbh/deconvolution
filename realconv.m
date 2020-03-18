function y = realconv(a,b)
if length(a) != length(b)
	printf 'Vectors are not of same length.';
else
	y = real(ifft(fft(a).*fft(b)));
end
