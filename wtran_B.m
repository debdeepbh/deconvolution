function [w] = wtran_B(z, type, p) 
% Copmute the wavelet coefficients using the basis matrix

% get the basis matrix
N = length(z);
B = getbasismat(type, p, N);

% empty wavelet coeff
w = [];

% Compute the wavelet transform using the basis matrix
for j=1:p+1
	Bsis = B(j,:);	% j-th row of basis matrix

	%%% Modification of the index for the coarsest level
	%%% will use indexnow as a replacement of j
	indexnow = j;
	if j == p+1
		indexnow = p;
	end

	% prepare beta(l), the l-th wavelet coefficient
	beta = zeros(1,N/(2^indexnow));

	for l=0:(N/(2^indexnow) - 1)	% range of l, the localization index for current level
		% this is Psi_{j,l}, the (j,l)-th basis element
		Psi = shift(Bsis,(2^indexnow)*l);

		% dot product with the Fourier-based deconvolution, the wavelet coefficient
		% matlab dot product, non-commutative, u.v = sum{conj(u),v}
		beta(l+1) = dot(Psi, z);
	end

	w = [w beta];
end




