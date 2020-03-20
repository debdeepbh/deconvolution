function [sigmal] = getleakedsd(B, L, r, noiseax, scaling) 
% computes the leaked noise s.d., given impulse responses and the 
% wavelet basis matrix B


[p_plus_one, N] = size(B);


for j=1:p_plus_one
	Bsis = B(j,:);	% j-th row of basis matrix

	to_mult = L ./ ((L + r * scaling(j)).^2);

	dotprod_val = sum( to_mult .* (abs(fft(Bsis)).^2));	% sum over the freq

	sigmal(j) = sqrt( noiseax^2/N * dotprod_val);
end
